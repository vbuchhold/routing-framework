#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <vector>

#include "Algorithms/Dijkstra/Dijkstra.h"
#include "DataStructures/Containers/BitVector.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"

// This class implements a closest-POI query in a road network, based on a simple Dijkstra search
// that stops when the k-th POI vertex is reached. It works in two phases. Given a set of POIs, the
// caller must first build a POI index by invoking the method buildPoiIndexFor(). The index is then
// used to run closest-POI queries by invoking the method findClosestPois().
template <typename GraphT, typename WeightT>
class DijkstraClosestPoiQuery {
 public:
  // Precomputed auxiliary information to accelerate POI queries.
  struct PoiIndex {
    friend class DijkstraClosestPoiQuery;

   public:
    // Returns the space (in bytes) consumed by this POI index.
    constexpr int spaceConsumption() const noexcept {
      return 0;
    }

   private:
    // Builds the POI index for the specified set of POI vertices.
    PoiIndex(const std::vector<int32_t>& pointsOfInterest, const int numVertices)
        : pointsOfInterest(pointsOfInterest), vertexContainsPoi(numVertices) {
      assert(!pointsOfInterest.empty());
      assert(std::is_sorted(pointsOfInterest.begin(), pointsOfInterest.end()));
      assert(numVertices > pointsOfInterest.back());
      for (const auto poi : pointsOfInterest)
        vertexContainsPoi[poi] = true;
    }

    const std::vector<int32_t>& pointsOfInterest; // The POI vertices in increasing order of ID.
    BitVector vertexContainsPoi;                  // Indicates for each vertex whether it is a POI.
  };

  // A point of interest returned by a POI query.
  struct Poi {
    int vertex; // The vertex that contains this POI.
    int dist;   // The shortest-path distance from the source to this POI.
  };

  // Creates an instance of a Dijkstra-based closest-POI query.
  DijkstraClosestPoiQuery(const GraphT& graph)
      : graph(graph),
        dijkstra(graph, {vertexContainsPoi, numRemainingPois, closestPois}),
        vertexContainsPoi(nullptr) {}

  // Builds the POI index for the specified set of POI vertices.
  PoiIndex buildPoiIndexFor(const std::vector<int32_t>& pointsOfInterest) const {
    return {pointsOfInterest, graph.numVertices()};
  }

  // Returns the k closest POI vertices to s.
  const std::vector<Poi>& findClosestPois(const int s, const PoiIndex& idx, const int k = 1) {
    assert(idx.vertexContainsPoi.size() == graph.numVertices());
    vertexContainsPoi = &idx.vertexContainsPoi;
    numRemainingPois = k;
    closestPois.clear();
    dijkstra.run(s);
    return closestPois;
  }

 private:
  // The stopping criterion for the Dijkstra search. We can stop when the k-th POI is reached.
  struct StopWhenClosestPoisReached {
    // Creates a stopping criterion for the Dijkstra search.
    StopWhenClosestPoisReached(
        const BitVector*& vertexContainsPoi, int& numRemainingPois,
        std::vector<Poi>& closestPois) noexcept
        : vertexContainsPoi(vertexContainsPoi),
          numRemainingPois(numRemainingPois),
          closestPois(closestPois) {}

    // Returns true if the search can be stopped at v.
    template <typename DistLabelT, typename DistLabelContainerT>
    bool operator()(const int v, DistLabelT& distToV, const DistLabelContainerT& /*distLabels*/) {
      assert(vertexContainsPoi != nullptr);
      if ((*vertexContainsPoi)[v]) {
        closestPois.push_back({v, distToV[0]});
        return --numRemainingPois == 0;
      }
      return false;
    }

    const BitVector*& vertexContainsPoi; // Indicates for each vertex whether it contains a POI.
    int& numRemainingPois;               // The number of POI vertices that are yet to be reached.
    std::vector<Poi>& closestPois;       // The closest POI vertices encountered so far.
  };

  using LabelSet = BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>;
  using DijkstraSearch = Dijkstra<GraphT, WeightT, LabelSet, StopWhenClosestPoisReached>;

  const GraphT& graph;                // The graph we work on.
  DijkstraSearch dijkstra;            // A standard Dijkstra search.
  const BitVector* vertexContainsPoi; // Indicates for each vertex whether it contains a POI.
  int numRemainingPois;               // The number of POI vertices that are yet to be reached.
  std::vector<Poi> closestPois;       // The closest POI vertices encountered so far.
};
