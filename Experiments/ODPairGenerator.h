#pragma once

#include <cassert>
#include <random>
#include <utility>
#include <vector>

#include "Algorithms/Dijkstra/Dijkstra.h"
#include "DataStructures/Containers/BitVector.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "DataStructures/Utilities/OriginDestination.h"
#include "Tools/Constants.h"
#include "Tools/ContainerHelpers.h"
#include "Tools/OpenMP.h"

// A facility generating O-D pairs (or queries) for the experimental evaluation of shortest-path
// algorithms. It supports choosing O and D uniformly at random, and more advanced methodologies
// such as choosing D by Dijkstra rank. Origin and destination vertices can be restricted to a
// smaller study area.
template <typename GraphT, typename WeightT>
class ODPairGenerator {
 public:
  // Constructs an OD-pair generator with the specified graph.
  ODPairGenerator(const GraphT& graph, const int seed = 0)
      : ODPairGenerator(graph, BitVector(graph.numVertices(), true), seed) {}

  // Constructs an OD-pair generator with the specified graph and study area.
  ODPairGenerator(const GraphT& graph, const BitVector& isVertexInStudyArea, const int seed = 0)
      : rand(seed + omp_get_thread_num() + 1),
        distribution(0, isVertexInStudyArea.cardinality() - 1),
        dijkstra(graph),
        isVertexInStudyArea(isVertexInStudyArea),
        verticesInStudyArea(isVertexInStudyArea.cardinality()) {
    assert(graph.numVertices() == isVertexInStudyArea.size());
    for (auto i = 0, u = isVertexInStudyArea.firstSetBit(); i < verticesInStudyArea.size(); ++i) {
      verticesInStudyArea[i] = u;
      u = isVertexInStudyArea.nextSetBit(u);
    }
  }

  // Returns an OD pair with O and D picked uniformly at random.
  OriginDestination getRandomODPair() {
    return {verticesInStudyArea[distribution(rand)], verticesInStudyArea[distribution(rand)]};
  }

  // Returns a random OD pair, where D has the specified Dijkstra rank from O.
  std::pair<OriginDestination, int> getRandomODPairChosenByRank(const int rank) {
    return getRandomODPairChosenByRank(rank, verticesInStudyArea[distribution(rand)]);
  }

  // Returns a random OD pair, where D has the specified Dijkstra rank from O.
  std::pair<OriginDestination, int> getRandomODPairChosenByRank(const int rank, const int src) {
    assert(rank >= 0);
    assert(contains(verticesInStudyArea.begin(), verticesInStudyArea.end(), src));
    auto dst = src;
    auto actualRank = -1;
    dijkstra.init({{src}});
    for (; !dijkstra.queue.empty() && actualRank < rank; actualRank += isVertexInStudyArea[dst])
      dst = dijkstra.settleNextVertex();
    return {{src, dst}, actualRank};
  }

  // Returns a random OD pair, where the distance between O and D is dist.
  std::pair<OriginDestination, int> getRandomODPairChosenByDistance(const int dist) {
    return getRandomODPairChosenByDistance(dist, verticesInStudyArea[distribution(rand)]);
  }

  // Returns a random OD pair, where the distance between O and D is dist.
  std::pair<OriginDestination, int> getRandomODPairChosenByDistance(const int dist, const int src) {
    assert(dist >= 0);
    assert(contains(verticesInStudyArea.begin(), verticesInStudyArea.end(), src));
    auto dst = src;
    dijkstra.init({{src}});
    while (!dijkstra.queue.empty() &&
           (dijkstra.getDistance(dst) < dist || !isVertexInStudyArea[dst]))
      dst = dijkstra.settleNextVertex();
    return {{src, dst}, dijkstra.getDistance(dst)};
  }

  // Returns the Dijkstra rank for the specified destination with respect to the given origin.
  int getDijkstraRankFor(const OriginDestination& od) {
    assert(contains(verticesInStudyArea.begin(), verticesInStudyArea.end(), od.origin));
    assert(contains(verticesInStudyArea.begin(), verticesInStudyArea.end(), od.destination));
    auto u = INVALID_VERTEX;
    auto rank = -1;
    dijkstra.init({{od.origin}});
    while (u != od.destination) {
      if (dijkstra.queue.empty())
        return INFTY;
      u = dijkstra.settleNextVertex();
      rank += isVertexInStudyArea[u];
    }
    return rank;
  }

 private:
  using Dijkstra = StandardDijkstra<GraphT, WeightT, BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>>;

  std::minstd_rand rand;                        // A Lehmer random number generator.
  std::uniform_int_distribution<> distribution; // A functor returning uniform random indices.
  Dijkstra dijkstra;                            // Dijkstra's shortest-path algorithm.

  const BitVector isVertexInStudyArea;  // Indicates whether a vertex is in the study area.
  std::vector<int> verticesInStudyArea; // All vertices in study area, i.e., all possible sources.
};
