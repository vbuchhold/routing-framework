#pragma once

#include <array>
#include <cassert>
#include <vector>

#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Labels/Containers/ParentLabelContainer.h"
#include "DataStructures/Labels/Containers/StampedDistanceLabelContainer.h"
#include "DataStructures/Queues/AddressableKHeap.h"
#include "Tools/Constants.h"

namespace dijkstra {

// A dummy pruning criterion for Dijkstra's algorithm that does no pruning at all.
struct NoPruningCriterion {
  // Returns always false.
  template <typename DistanceLabelT, typename DistanceLabelContainerT>
  bool operator()(const int, const DistanceLabelT&, DistanceLabelContainerT&) const {
    return false;
  }
};

}

// Implementation of Dijkstra's shortest-path algorithm. Depending on the used label set, it keeps
// track of parent vertices and/or edges, and computes multiple shortest paths simultaneously,
// possibly using SSE or AVX instructions. The algorithm can be used with different distance label
// containers and priority queues.
template <
    typename GraphT, typename WeightT, template <typename> class DistanceLabelContainerT,
    typename LabelSetT, typename QueueT, typename PruningCriterionT = dijkstra::NoPruningCriterion>
class Dijkstra {
  // Some classes are allowed to execute a Dijkstra search step by step.
  template <typename, template <typename> class>
  friend class BiDijkstra;
  template <typename, typename, typename>
  friend class ODPairGenerator;

 private:
  using Graph    = GraphT;                    // The graph type on which we compute shortest paths.
  using LabelSet = LabelSetT;                 // The distance and parent label type we use.
  using Queue    = QueueT;                    // The priority queue type.
  using PruningCriterion = PruningCriterionT; // The criterion applied to prune the search.

  using LabelMask = typename LabelSet::LabelMask;         // Marks subset of components in a label.
  using DistanceLabel = typename LabelSet::DistanceLabel; // The distance label of a vertex.
  using ParentLabel = typename LabelSet::ParentLabel;     // The parent label of a vertex.

  static constexpr int K = LabelSet::K; // The number of simultaneous shortest-path computations.

 public:
  // Constructs a Dijkstra instance.
  Dijkstra(const Graph& graph, const PruningCriterionT& pruneSearch = {})
      : graph(graph),
        distanceLabels(graph.numVertices()),
        parent(graph),
        queue(graph.numVertices()),
        pruneSearch(pruneSearch) {}

  // Ensures that the internal data structures fit for the size of the graph.
  void resize() {
    distanceLabels.resize(graph.numVertices());
    parent.resize();
    queue.resize(graph.numVertices());
  }

  // Runs a Dijkstra search from s to t. If t is omitted, runs a one-to-all search from s.
  void run(const int s, const int t = INVALID_VERTEX) {
    std::array<int, K> sources;
    std::fill(sources.begin(), sources.end(), s);
    init(sources);
    while (!queue.empty()) {
      const int settled = settleNextVertex();
      if (settled == t)
        break;
    }
  }

  // Runs a Dijkstra search that computes multiple shortest paths simultaneously.
  void run(const std::array<int, K>& sources, const std::array<int, K>& targets) {
    init(sources);
    while (!queue.empty()) {
      // Stop the search as soon as Q.minKey >= d_i(t_i) for all i = 1, ..., k.
      bool stop = true;
      for (int i = 0; i < K; ++i)
        stop &= queue.minKey() >= distanceLabels[targets[i]][i];
      if (stop)
        break;
      settleNextVertex();
    }
  }

  // Returns the shortest-path distance from the i-th source to t.
  int getDistance(const int t, const int i = 0) {
    return distanceLabels[t][i];
  }

  // Returns the vertices along the shortest path from the i-th source to t in reverse order.
  std::vector<int> getReversePath(int t, const int i = 0) {
    assert(distanceLabels[t][i] != INFTY);
    return parent.getReversePath(t, i);
  }

  // Returns the edges along the shortest path from the i-th source to t in reverse order.
  std::vector<int> getReverseEdgePath(int t, const int i = 0) {
    assert(distanceLabels[t][i] != INFTY);
    return parent.getReverseEdgePath(t, i);
  }

 private:
  // Resets the distance labels and inserts all k simultaneous sources into the queue.
  void init(const std::array<int, K>& sources) {
    distanceLabels.init();
    queue.clear();
    for (int i = 0; i < K; ++i) {
      const int s = sources[i];
      distanceLabels[s][i] = 0;
      parent.setVertex(s, INVALID_VERTEX, true);
      parent.setEdge(s, INVALID_EDGE, true);
      if (!queue.contains(s))
        queue.insert(s, 0);
    }
  }

  // Removes the next vertex from the queue, relaxes its outgoing edges and returns its ID.
  int settleNextVertex() {
    int u, key;
    queue.deleteMin(u, key);
    const DistanceLabel& distToU = distanceLabels[u];

    // Check if the search can be pruned at u.
    if (pruneSearch(u, distToU, distanceLabels))
      return u;

    // Relax all edges out of u.
    FORALL_INCIDENT_EDGES(graph, u, e) {
      const int v = graph.edgeHead(e);
      DistanceLabel& distToV = distanceLabels[v];
      const DistanceLabel tentativeDist = distToU + graph.template get<WeightT>(e);
      const LabelMask mask = tentativeDist < distToV;
      if (mask) {
        distToV.min(tentativeDist);
        parent.setVertex(v, u, mask);
        parent.setEdge(v, e, mask);
        if (queue.contains(v))
          queue.decreaseKey(v, distToV.getKey());
        else
          queue.insert(v, distToV.getKey());
      }
    }
    return u;
  }

  using DistanceLabelCont = DistanceLabelContainerT<DistanceLabel>;
  using ParentLabelCont   = ParentLabelContainer<Graph, LabelSet>;

  const Graph& graph;               // The graph on which we compute shortest paths.
  DistanceLabelCont distanceLabels; // The distance labels of the vertices.
  ParentLabelCont parent;           // The parent information for each vertex.
  Queue queue;                      // The priority queue of unsettled vertices.
  PruningCriterionT pruneSearch;    // A criterion applied to prune the search.
};

// An alias template for a standard Dijkstra search.
template <typename GraphT, typename WeightT, typename LabelSetT>
using StandardDijkstra =
    Dijkstra<GraphT, WeightT, StampedDistanceLabelContainer, LabelSetT, AddressableQuadheap>;
