#pragma once

#include <cassert>
#include <random>

#include "Algorithms/Dijkstra/Dijkstra.h"
#include "DataStructures/Graph/Attributes/TravelTimeAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "Tools/Constants.h"
#include "Tools/OpenMP.h"

// A Dijkstra-based approach to choose the closest opportunity with high fitness. In the radiation
// model, each opportunity has a fitness value drawn from some distribution, and each traveler has
// a threshold drawn from the same distribution. Then the traveler chooses the closest opportunity
// with a fitness value higher than their threshold.
template <typename GraphT>
class DijkstraOpportunityChooser {
 public:
  // Constructs an opportunity chooser based on Dijkstra's algorithm.
  explicit DijkstraOpportunityChooser(const GraphT& graph, const int seed)
      : graph(graph), numOpportunities(0), rand(seed + omp_get_thread_num() + 1), dijkstra(graph) {
    FORALL_VERTICES(graph, u)
      numOpportunities += graph.numOpportunities(u);
    assert(numOpportunities > 0);
    assert(seed >= 0);
  }

  // Returns the vertex with the closest opportunity with sufficiently high fitness.
  int findClosestOpportunityWithHighFitness(const int src, const int numFitOpportunities) {
    assert(numFitOpportunities > 0);
    std::geometric_distribution<> dist(1.0 * numFitOpportunities / numOpportunities);
    int u = INVALID_VERTEX, numInterveningOpportunities = numOpportunities;
    while (numInterveningOpportunities >= numOpportunities)
      numInterveningOpportunities = dist(rand);
    dijkstra.init({{src}});
    while (numInterveningOpportunities >= 0) {
      assert(!dijkstra.queue.empty());
      u = dijkstra.settleNextVertex();
      numInterveningOpportunities -= graph.numOpportunities(u);
    }
    return u;
  }

  // Returns the cumulative number of opportunities between source and target.
  int computeInterveningOpportunities(const int src, const int dst) noexcept {
    int u = INVALID_VERTEX;
    int numInterveningOpportunities = -graph.numOpportunities(src) - graph.numOpportunities(dst);
    dijkstra.init({{src}});
    while (u != dst) {
      assert(!dijkstra.queue.empty());
      u = dijkstra.settleNextVertex();
      numInterveningOpportunities += graph.numOpportunities(u);
    }
    assert(numInterveningOpportunities >= 0);
    return numInterveningOpportunities;
  }

 private:
  using LabelSet = BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>;
  using Dijkstra = StandardDijkstra<GraphT, TravelTimeAttribute, LabelSet>;

  const GraphT& graph;   // The network we work on.
  int numOpportunities;  // The total number of opportunities in the network.
  std::minstd_rand rand; // A Lehmer random number generator.
  Dijkstra dijkstra;     // Dijkstra's shortest-path algorithm.
};
