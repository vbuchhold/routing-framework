#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <random>
#include <vector>

#include <stocc/stocc.h>

#include "Algorithms/CCH/CCH.h"
#include "Algorithms/CCH/EliminationTreeQuery.h"
#include "Algorithms/CH/CH.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "Tools/Constants.h"
#include "Tools/OpenMP.h"

// A CCH-based approach to choose the closest opportunity with high fitness. In the radiation
// model, each opportunity has a fitness value drawn from some distribution, and each traveler has
// a threshold drawn from the same distribution. Then the traveler chooses the closest opportunity
// with a fitness value higher than their threshold.
template <typename GraphT>
class CCHOpportunityChooser {
 public:
  // Creates an opportunity chooser based on the specified CCH (and the corresponding minimum CH).
  CCHOpportunityChooser(const GraphT& graph, const int seed, const CCH& cch, const CH& minCH)
      : urbg(seed + omp_get_thread_num() + 1),
        nrng(seed + omp_get_thread_num()),
        cch(cch),
        elimTreeQuery(minCH, cch.getEliminationTree()),
        numLocationsAmongFirst(graph.numVertices() + 1),
        numOpportunitiesAmongFirst(graph.numVertices() + 1) {
    assert(graph.numVertices() == cch.getUpwardGraph().numVertices());
    assert(graph.numVertices() == minCH.upwardGraph().numVertices());
    assert(seed >= 0);

    auto numLocations = 0;
    for (auto i = 0; i < graph.numVertices(); ++i)
      numLocations += graph.numOpportunities(i) > 0;
    locations.assign(numLocations, 0);
    numOpportunitiesPerLocation.assign(numLocations, 0);

    for (auto i = 0, j = 0; i < graph.numVertices(); ++i) {
      const auto numOpportunities = graph.numOpportunities(cch.getContractionOrder()[i]);
      if (numOpportunities > 0) {
        locations[j] = i;
        numOpportunitiesPerLocation[j] = numOpportunities;
        ++j;
      }
      numLocationsAmongFirst[i + 1] = j;
      numOpportunitiesAmongFirst[i + 1] = numOpportunitiesAmongFirst[i] + numOpportunities;
    }

    opportunityCounts.assign(
        std::max(cch.getSeparatorDecomposition().tree.size(), locations.size()), 0);
  }

  // Returns the vertex with the closest opportunity with sufficiently high fitness.
  int findClosestOpportunityWithHighFitness(const int src, const int numFitOpportunities) {
    assert(numFitOpportunities > 0);
    const auto& decomp = cch.getSeparatorDecomposition();
    const auto s = cch.getRanks()[src];
    elimTreeQuery.pinForwardSearch(s);
    recursionStack.push_back({0, 0, 0, numFitOpportunities});
    distToClosestOpportunity = INFTY;

    while (!recursionStack.empty()) {
      const auto parent = recursionStack.back();
      recursionStack.pop_back();

      if (parent.dist >= distToClosestOpportunity)
        continue;

      const auto firstLoc = numLocationsAmongFirst[parent.firstVertex];
      const auto lastLoc = numLocationsAmongFirst[decomp.lastSeparatorVertex(parent.node)];
      if (parent.numFitOpportunities * int64_t{lastLoc - firstLoc} <= RECURSION_THRESHOLD) {
        const auto firstOp = numOpportunitiesAmongFirst[parent.firstVertex];
        const auto lastOp = numOpportunitiesAmongFirst[decomp.lastSeparatorVertex(parent.node)];
        handleBaseCase(firstLoc, lastLoc, parent.numFitOpportunities, lastOp - firstOp);
        continue;
      }

      // Retrieve the number of opportunities in each child region.
      auto child = decomp.leftChild(parent.node), firstVertex = parent.firstVertex, i = 0;
      while (child != 0) {
        const auto firstOpInChild = numOpportunitiesAmongFirst[firstVertex];
        const auto lastOpInChild = numOpportunitiesAmongFirst[decomp.lastSeparatorVertex(child)];
        opportunityCounts[i++] = lastOpInChild - firstOpInChild;
        firstVertex = decomp.lastSeparatorVertex(child);
        child = decomp.rightSibling(child);
      }
      assert(firstVertex == decomp.firstSeparatorVertex(parent.node));

      // Retrieve the number of opportunities in the separator.
      const auto firstOpInSep = numOpportunitiesAmongFirst[firstVertex];
      const auto lastOpInSep = numOpportunitiesAmongFirst[decomp.lastSeparatorVertex(parent.node)];
      opportunityCounts[i++] = lastOpInSep - firstOpInSep;

      // Draw the number of fit opportunities in each child region and in the separator.
      const auto numDraws = parent.numFitOpportunities;
      nrng.MultiHypergeometric(opportunityCounts.data(), opportunityCounts.data(), numDraws, i);

      // Push each child region with a nonzero number of fit opportunities onto the recursion stack.
      const auto stackSize = recursionStack.size();
      child = decomp.leftChild(parent.node), firstVertex = parent.firstVertex, i = 0;
      while (child != 0) {
        if (opportunityCounts[i] > 0) {
          const auto childContainsSrc = firstVertex <= s && s < decomp.lastSeparatorVertex(child);
          const auto dist = childContainsSrc ? 0 : computeDistToSubgraph(child);
          recursionStack.push_back({child, dist, firstVertex, opportunityCounts[i]});
        }
        ++i;
        firstVertex = decomp.lastSeparatorVertex(child);
        child = decomp.rightSibling(child);
      }
      assert(firstVertex == decomp.firstSeparatorVertex(parent.node));
      std::sort(recursionStack.rbegin(), recursionStack.rend() - stackSize);

      // Sample fit opportunities from the separator (if it contains any fit opportunities).
      if (opportunityCounts[i] > 0) {
        const auto firstLocInSep = numLocationsAmongFirst[firstVertex];
        handleBaseCase(firstLocInSep, lastLoc, opportunityCounts[i], lastOpInSep - firstOpInSep);
      }
    }

    assert(distToClosestOpportunity != INFTY);
    return cch.getContractionOrder()[closestOpportunity];
  }

 private:
  // A subgraph of the hierarchy that has not been examined yet.
  struct ActiveSubgraph {
    // Returns true if this subgraph is closer to the source than the specified subgraph.
    bool operator<(const ActiveSubgraph& rhs) const noexcept {
      return dist < rhs.dist;
    }

    int32_t node;                // The separator decomposition node corresponding to this subgraph.
    int32_t dist;                // A lower bound on the distance to any vertex in this subgraph.
    int32_t firstVertex;         // The first vertex in this subgraph.
    int32_t numFitOpportunities; // The number of sufficiently fit opportunities in this subgraph.
  };

  static constexpr int RECURSION_THRESHOLD = 4096; // Used to stop the recursion during queries.

  // Samples numFitOpportunities opportunities from the specified range of locations and checks for
  // each whether it improves the closest fit opportunity so far encountered.
  void handleBaseCase(
      const int firstLoc, const int lastLoc,
      const int numFitOpportunities, const int numOpportunities) {
    assert(0 <= firstLoc); assert(firstLoc <= lastLoc); assert(lastLoc <= locations.size());
    assert(numFitOpportunities <= numOpportunities);
    const auto beginning = numOpportunitiesPerLocation.begin();
    std::copy(beginning + firstLoc, beginning + lastLoc, opportunityCounts.begin());
    for (auto i = 1; i <= numFitOpportunities; ++i) {
      const auto rank = std::uniform_int_distribution<>(0, numOpportunities - i)(urbg);
      auto j = 0;
      auto partialSum = opportunityCounts.front();
      while (partialSum <= rank) partialSum += opportunityCounts[++j];
      assert(j < lastLoc - firstLoc);
      assert(partialSum <= numOpportunities);
      assert(opportunityCounts[j] > 0);
      --opportunityCounts[j];
      const auto loc = locations[firstLoc + j];
      elimTreeQuery.runReverseSearch(loc);
      const auto dist = elimTreeQuery.getDistance();
      if (dist < distToClosestOpportunity) {
        closestOpportunity = loc;
        distToClosestOpportunity = dist;
      }
    }
  }

  // Returns a lower bound on the distance from the source to any vertex in the specified subgraph.
  int computeDistToSubgraph(const int node) {
    const auto& upwardGraph = cch.getUpwardGraph();
    const auto highestVertex = cch.getSeparatorDecomposition().lastSeparatorVertex(node) - 1;
    const auto firstBoundaryVertex = &upwardGraph.edgeHead(upwardGraph.firstEdge(highestVertex));
    const auto lastBoundaryVertex = &upwardGraph.edgeHead(upwardGraph.lastEdge(highestVertex));
    elimTreeQuery.runReverseSearch(firstBoundaryVertex, lastBoundaryVertex);
    return elimTreeQuery.getDistance();
  }

  using ElimTreeQuery = EliminationTreeQuery<BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>>;

  std::minstd_rand urbg; // A uniform random bit generator.
  StochasticLib1 nrng;   // A nonuniform random number generator.

  const CCH& cch;              // The metric-independent CCH.
  ElimTreeQuery elimTreeQuery; // A standard elimination tree query.

  std::vector<int32_t> locations;                   // Vertices with nonzero no. of opportunities.
  std::vector<int32_t> numOpportunitiesPerLocation; // No. of opportunities for each location.

  std::vector<int32_t> numLocationsAmongFirst;     // No. of locations among first i vertices.
  std::vector<int32_t> numOpportunitiesAmongFirst; // No. of opportunities among first i vertices.

  std::vector<ActiveSubgraph> recursionStack; // The stack of active subgraphs.
  std::vector<int32_t> opportunityCounts;     // Temporary storage holding opportunity counts.
  int closestOpportunity;                     // The closest fit opportunity so far encountered.
  int distToClosestOpportunity;               // The distance to the closest fit opportunity.
};
