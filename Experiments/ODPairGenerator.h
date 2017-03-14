#pragma once

#include <cassert>
#include <random>

#include "Algorithms/Dijkstra/Dijkstra.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Utilities/OriginDestination.h"
#include "Tools/Constants.h"

// A facility generating O-D pairs (or queries) for the experimental evaluation of shortest-path
// algorithms. It supports choosing O and D uniformly at random, and more advanced methodologies
// such as choosing D by Dijkstra rank.
template <
    typename GraphT, typename WeightT, typename RandomNumberGeneratorT = std::default_random_engine>
class ODPairGenerator {
 public:
  // Constructs an O-D pair generator for the specified graph.
  ODPairGenerator(const GraphT& graph, RandomNumberGeneratorT& rand)
      : dijkstra(graph),
        rand(rand),
        dist(0, graph.numVertices() - 1) {}

  // Returns an O-D pair with O and D picked uniformly at random.
  OriginDestination getRandomODPair() {
    return {dist(rand), dist(rand)};
  }

  // Returns a random O-D pair, where D has the specified Dijkstra rank from O.
  OriginDestination getRandomODPairChosenByDijkstraRank(const int rank) {
    assert(rank >= 0); assert(rank <= dist.b());
    const int o = dist(rand);
    int d = o;
    dijkstra.init({o});
    for (int i = 0; i <= rank; ++i) {
      if (dijkstra.queue.empty())
        assert(false);
      d = dijkstra.settleNextVertex();
    }
    return {o, d};
  }

  // Returns a random O-D pair, where the distance between O and D is distance.
  OriginDestination getRandomODPairChosenByDistance(const int distance) {
    assert(distance >= 0);
    const int o = dist(rand);
    int d = o;
    dijkstra.init({o});
    while (!dijkstra.queue.empty() && dijkstra.getDistance(d) < distance)
      d = dijkstra.settleNextVertex();
    return {o, d};
  }

 private:
  using LabelSet = BasicLabelSet<1, ParentInfo::NO_PARENT_INFO>;
  using Dijkstra = StandardDijkstra<GraphT, WeightT, LabelSet>;

  Dijkstra dijkstra;                    // The Dijkstra search used to compute Dijkstra ranks.
  RandomNumberGeneratorT& rand;         // The random number engine providing randomness.
  std::uniform_int_distribution<> dist; // A functor returning random vertex IDs.
};
