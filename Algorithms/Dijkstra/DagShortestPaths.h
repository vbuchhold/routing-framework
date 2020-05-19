#pragma once

#include <cassert>
#include <cstdint>
#include <vector>

#include "Algorithms/Dijkstra/Dijkstra.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Labels/Containers/ParentLabelContainer.h"
#include "DataStructures/Labels/Containers/StampedDistanceLabelContainer.h"
#include "DataStructures/Queues/AddressableKHeap.h"
#include "Tools/Constants.h"

// Implementation of a shortest-path search on a directed acyclic graph. The vertices must be
// numbered in topological order. The search works similarly to Dijkstra's algorithm, but processes
// vertices in topological order rather than in increasing order of distance, and thus needs no
// decrease-key operations on the priority queue. Depending on the used label set, it keeps parent
// vertices and/or edges. The search can be used with different distance label containers and
// priority queues. Moreover, the caller can provide an own pruning criterion.
template <
    typename GraphT, typename WeightT, typename LabelSetT,
    typename PruningCriterionT = dij::NoCriterion,
    template <typename> class DistanceLabelContainerT = StampedDistanceLabelContainer,
    typename QueueT = AddressableQuadheap>
class DagShortestPaths {
 public:
  // Constructs a shortest-path search instance for the specified directed acyclic graph.
  explicit DagShortestPaths(const GraphT& graph, PruningCriterionT pruneSearch = {})
      : graph(graph),
        distanceLabels(graph.numVertices()),
        parent(graph),
        queue(graph.numVertices()),
        pruneSearch(pruneSearch) {}

  // Runs a shortest-path search from s.
  void run(const int s) {
    runWithOffset(s, 0);
  }

  // Runs a shortest-path search from s to t.
  void run(const int s, const int t) {
    init(s);
    while (!queue.empty()) {
      if (queue.minId() == t)
        break;
      settleNextVertex();
    }
  }

  // Runs a shortest-path search from s, with the distance of s initialized to the given offset.
  void runWithOffset(const int s, const int offset) {
    init(s, offset);
    while (!queue.empty())
      settleNextVertex();
  }

  // Returns the shortest-path distance to t.
  int getDistance(const int t) {
    return distanceLabels[t][0];
  }

  // Returns the parent vertex of v on the shortest path to v.
  int getParentVertex(const int v) {
    assert(distanceLabels[v][0] != INFTY);
    return parent.getVertex(v);
  }

  // Returns the parent edge of v on the shortest path to v.
  int getParentEdge(const int v) {
    assert(distanceLabels[v][0] != INFTY);
    return parent.getEdge(v);
  }

  // Returns the vertices on the shortest path to t in reverse order.
  const std::vector<int32_t>& getReversePath(const int t) {
    assert(distanceLabels[t][0] != INFTY);
    return parent.getReversePath(t);
  }

  // Returns the edges on the shortest path to t in reverse order.
  const std::vector<int32_t>& getReverseEdgePath(const int t) {
    assert(distanceLabels[t][0] != INFTY);
    return parent.getReverseEdgePath(t);
  }

 private:
  // Resets the distance labels and inserts the source into the queue.
  void init(const int s, const int offset = 0) {
    distanceLabels.init();
    queue.clear();
    distanceLabels[s] = offset;
    parent.setVertex(s, s, true);
    parent.setEdge(s, INVALID_EDGE, true);
    queue.insert(s, s);
  }

  // Removes the next vertex from the queue, relaxes its outgoing edges, and returns its ID.
  int settleNextVertex() {
    int v, key;
    queue.deleteMin(v, key);
    auto& distToV = distanceLabels[v];

    // Check whether the search can be pruned at v.
    if (pruneSearch(v, distToV, distanceLabels))
      return v;

    // Relax all edges out of v.
    FORALL_INCIDENT_EDGES(graph, v, e) {
      const auto w = graph.edgeHead(e);
      auto& distToW = distanceLabels[w];
      const auto distViaV = distToV + graph.template get<WeightT>(e);
      const auto mask = distViaV < distToW;
      if (mask) {
        distToW.min(distViaV);
        parent.setVertex(w, v, mask);
        parent.setEdge(w, e, mask);
        if (!queue.contains(w))
          queue.insert(w, w);
      }
    }
    return v;
  }

  using DistanceLabelCont = DistanceLabelContainerT<typename LabelSetT::DistanceLabel>;
  using ParentLabelCont = ParentLabelContainer<GraphT, LabelSetT>;

  const GraphT& graph;              // The graph (DAG) on which we compute shortest paths.
  DistanceLabelCont distanceLabels; // The distance labels of the vertices.
  ParentLabelCont parent;           // The parent information for each vertex.
  QueueT queue;                     // The priority queue of unsettled vertices.
  PruningCriterionT pruneSearch;    // The criterion used to prune the search.
};
