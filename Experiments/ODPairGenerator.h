#pragma once

#include <cassert>
#include <random>
#include <utility>

#include "Algorithms/Dijkstra/Dijkstra.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "DataStructures/Utilities/OriginDestination.h"
#include "Tools/Constants.h"
#include "Tools/OpenMP.h"

// A facility generating O-D pairs (or queries) for the experimental evaluation of shortest-path
// algorithms. It supports choosing O and D uniformly at random, and more advanced methodologies
// such as choosing D by Dijkstra rank.
template <typename GraphT, typename WeightT>
class ODPairGenerator {
 public:
  // Constructs an OD-pair generator for the specified graph.
  ODPairGenerator(const GraphT& graph, const int seed = 0)
      : rand(seed + omp_get_thread_num() + 1),
        distribution(0, graph.numVertices() - 1),
        dijkstra(graph) {}

  // Returns an OD pair with O and D picked uniformly at random.
  OriginDestination getRandomODPair() {
    return {distribution(rand), distribution(rand)};
  }

  // Returns a random OD pair, where D has the specified Dijkstra rank from O.
  std::pair<OriginDestination, int> getRandomODPairChosenByRank(const int rank) {
    return getRandomODPairChosenByRank(rank, distribution(rand));
  }

  // Returns a random OD pair, where D has the specified Dijkstra rank from O.
  std::pair<OriginDestination, int> getRandomODPairChosenByRank(const int rank, const int src) {
    assert(rank >= 0);
    assert(src >= 0); assert(src <= distribution.b());
    auto dst = src;
    auto actualRank = -1;
    dijkstra.init({{src}});
    for (; !dijkstra.queue.empty() && actualRank < rank; ++actualRank)
      dst = dijkstra.settleNextVertex();
    return {{src, dst}, actualRank};
  }

  // Returns a random OD pair, where the distance between O and D is dist.
  std::pair<OriginDestination, int> getRandomODPairChosenByDistance(const int dist) {
    return getRandomODPairChosenByDistance(dist, distribution(rand));
  }

  // Returns a random OD pair, where the distance between O and D is dist.
  std::pair<OriginDestination, int> getRandomODPairChosenByDistance(const int dist, const int src) {
    assert(dist >= 0);
    assert(src >= 0); assert(src <= distribution.b());
    auto dst = src;
    dijkstra.init({{src}});
    while (!dijkstra.queue.empty() && dijkstra.getDistance(dst) < dist)
      dst = dijkstra.settleNextVertex();
    return {{src, dst}, dijkstra.getDistance(dst)};
  }

  // Returns the Dijkstra rank for the specified destination with respect to the given origin.
  int getDijkstraRankFor(const OriginDestination& od) {
    assert(od.origin >= 0); assert(od.origin <= distribution.b());
    assert(od.destination >= 0); assert(od.destination <= distribution.b());
    auto rank = 0;
    dijkstra.init({{od.origin}});
    while (dijkstra.settleNextVertex() != od.destination) {
      if (dijkstra.queue.empty())
        return INFTY;
      ++rank;
    }
    return rank;
  }

 private:
  using Dijkstra = StandardDijkstra<GraphT, WeightT, BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>>;

  std::minstd_rand rand;                        // A Lehmer random number generator.
  std::uniform_int_distribution<> distribution; // A functor returning uniform random vertex IDs.
  Dijkstra dijkstra;                            // Dijkstra's shortest-path algorithm.
};
