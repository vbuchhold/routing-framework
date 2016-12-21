#pragma once

#include <array>
#include <cassert>
#include <type_traits>
#include <vector>

#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Labels/Containers/StampedDistanceLabelContainer.h"
#include "DataStructures/Queues/Heap.h"
#include "Tools/Constants.h"

// Implementation of Dijkstra's shortest-path algorithm. Depending on the used label set, it keeps
// track of parent vertices and/or edges, and computes multiple shortest paths simultaneously,
// possibly using SSE or AVX instructions. The algorithm can be used with different distance label
// containers and priority queues.
template <
    typename GraphT, template <typename> class DistanceLabelContainerT, typename LabelSetT,
    typename QueueT, template <typename> class GetWeightT>
class Dijkstra {
  // Bidirectional Dijkstra is allowed to execute a Dijkstra search step by step.
  template <typename DijkstraT>
  friend class BiDijkstra;

 private:
  using Graph = GraphT;                                    // The graph type on which we operate.
  using LabelMask = typename LabelSetT::LabelMask;         // Marks subset of components in a label.
  using DistanceLabel = typename LabelSetT::DistanceLabel; // The distance label of a vertex.

  static constexpr int K = LabelSetT::K; // The number of simultaneous shortest-path computations.

 public:
  // Constructs a Dijkstra instance.
  Dijkstra(const Graph& graph)
      : graph(graph),
        distanceLabels(graph.numVertices()),
        parent(LabelSetT::KEEP_PARENT_VERTICES ? graph.numVertices() : 0),
        queue(graph.numVertices()) {}

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
  template <bool cond = LabelSetT::KEEP_PARENT_VERTICES, typename = std::enable_if_t<cond>>
  std::vector<int> getReversePath(int t, const int i = 0) {
    assert(distanceLabels[t][i] != INFTY);
    std::vector<int> path;
    while (parent[t].vertex(i) != INVALID_VERTEX) {
      assert(graph.containsEdge(parent[t].vertex(i), t));
      path.push_back(t);
      t = parent[t].vertex(i);
    }
    path.push_back(t);
    return path;
  }

  // Returns the edges along the shortest path from the i-th source to t in reverse order.
  template <bool cond = LabelSetT::KEEP_PARENT_EDGES, typename = std::enable_if_t<cond>>
  std::vector<int> getReverseEdgePath(int t, const int i = 0) {
    assert(distanceLabels[t][i] != INFTY);
    std::vector<int> path;
    while (parent[t].vertex(i) != INVALID_VERTEX) {
      assert(graph.containsEdge(parent[t].vertex(i), t));
      assert(graph.edgeHead(parent[t].edge(i)) == t);
      path.push_back(parent[t].edge(i));
      t = parent[t].vertex(i);
    }
    return path;
  }

 private:
  // Resets the distance labels and inserts all k simultaneous sources into the queue.
  void init(const std::array<int, K>& sources) {
    distanceLabels.init();
    queue.clear();
    for (int i = 0; i < K; ++i) {
      const int s = sources[i];
      distanceLabels[s][i] = 0;
      setParentVertex(s, INVALID_VERTEX, LabelMask(i));
      setParentEdge(s, INVALID_EDGE, LabelMask(i));
      if (!queue.contains(s))
        queue.insert(s, 0);
    }
  }

  // Removes the next vertex from the queue, relaxes its outgoing edges and returns its ID.
  int settleNextVertex() {
    int u, key;
    queue.deleteMin(u, key);
    const DistanceLabel& distToU = distanceLabels[u];

    // Relax all edges out of u.
    FORALL_INCIDENT_EDGES(graph, u, e) {
      const int v = graph.edgeHead(e);
      DistanceLabel& distToV = distanceLabels[v];
      const DistanceLabel tentativeDist = distToU + getWeight(graph, e);
      const LabelMask mask = tentativeDist < distToV;
      if (mask) {
        distToV.min(tentativeDist);
        setParentVertex(v, u, mask);
        setParentEdge(v, e, mask);
        if (queue.contains(v))
          queue.decreaseKey(v, distToV.getKey());
        else
          queue.insert(v, distToV.getKey());
      }
    }
    return u;
  }

  // Sets the parent vertex of v to u on all shortest paths specified by mask.

  template <bool cond = LabelSetT::KEEP_PARENT_VERTICES>
  std::enable_if_t<cond> setParentVertex(const int v, const int u, const LabelMask& mask) {
    assert(v >= 0); assert(v < parent.size());
    parent[v].setVertex(u, mask);
  }

  template <bool cond = LabelSetT::KEEP_PARENT_VERTICES>
  std::enable_if_t<!cond> setParentVertex(const int, const int, const LabelMask&) {}

  // Sets the parent edge of v to e on all shortest paths specified by mask.

  template <bool cond = LabelSetT::KEEP_PARENT_EDGES>
  std::enable_if_t<cond> setParentEdge(const int v, const int e, const LabelMask& mask) {
    assert(v >= 0); assert(v < parent.size());
    parent[v].setEdge(e, mask);
  }

  template <bool cond = LabelSetT::KEEP_PARENT_EDGES>
  std::enable_if_t<!cond> setParentEdge(const int, const int, const LabelMask&) {}

  using DistanceLabelContainer = DistanceLabelContainerT<typename LabelSetT::DistanceLabel>;
  using ParentLabelContainer = std::vector<typename LabelSetT::ParentLabel>;

  const Graph& graph;                    // The graph on which we compute shortest paths.
  DistanceLabelContainer distanceLabels; // The distance labels of the vertices.
  ParentLabelContainer parent;           // The parent information for each vertex.
  QueueT queue;                          // The priority queue of unsettled vertices.
  GetWeightT<Graph> getWeight;           // A functor returning the edge weight used for routing.
};

// An alias template for a standard Dijkstra search.
template <typename GraphT, typename LabelSetT, template <typename> class GetWeightT>
using StandardDijkstra =
    Dijkstra<GraphT, StampedDistanceLabelContainer, LabelSetT, QuadHeap, GetWeightT>;
