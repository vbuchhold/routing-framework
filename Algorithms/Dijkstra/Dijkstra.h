#pragma once

#include <cassert>
#include <type_traits>
#include <vector>

#include "DataStructures/DistanceLabels/StampedDistanceLabelContainer.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Queues/Heap.h"
#include "Tools/Constants.h"

// Implementation of Dijkstra's well-known shortest-path algorithm. When enabled by template
// arguments, it keeps track of parent vertices and/or edges. The algorithm can be used with
// different distance labels and priority queues.
template <
    typename GraphT, typename DistanceLabelContainerT, typename QueueT,
    template <typename> class GetWeightT, bool keepParentVertices, bool keepParentEdges>
class Dijkstra {
 public:
  // Constructs a Dijkstra instance.
  Dijkstra(const GraphT& graph)
      : graph(graph),
        distanceLabels(graph.numVertices()),
        parent(keepParentVertices || keepParentEdges ? graph.numVertices() : 0),
        queue(graph.numVertices()) {}

  // Runs a Dijkstra search from s to t. If t is omitted, runs a one-to-all query.
  void run(const int s, const int t = INVALID_VERTEX) {
    // Initialization.
    distanceLabels.init();
    distanceLabels[s] = 0;
    setParentVertex(s, INVALID_VERTEX);
    setParentEdge(s, INVALID_EDGE);
    queue.clear();
    queue.insert(s, 0);

    while (!queue.empty()) {
      // Settle u.
      int u, distToU;
      queue.deleteMin(u, distToU);
      if (u == t)
        break;

      // Relax all edges out of u.
      FORALL_INCIDENT_EDGES(graph, u, e) {
        const int v = graph.edgeHead(e);
        int& distToV = distanceLabels[v];
        if (distToU + getWeight(graph, e) < distToV) {
          distToV = distToU + getWeight(graph, e);
          setParentVertex(v, u);
          setParentEdge(v, e);
          if (queue.contains(v))
            queue.decreaseKey(v, distToV);
          else
            queue.insert(v, distToV);
        }
      }
    }
  }

  // Returns the shortest-path distance to t.
  int getDistance(const int t) {
    return distanceLabels[t];
  }

  // Returns the sequence of vertices along the shortest path to t in reverse order.
  template <bool cond = keepParentVertices || keepParentEdges, typename = std::enable_if_t<cond>>
  std::vector<int> getReversePath(int t) {
    assert(distanceLabels[t] != INFTY);
    std::vector<int> path;

    while (parent[t].vertex != INVALID_VERTEX) {
      assert(graph.containsEdge(parent[t].vertex, t));
      path.push_back(t);
      t = parent[t].vertex;
    }
    path.push_back(t);
    return path;
  }

  // Returns the sequence of edges along the shortest path to t in reverse order.
  template <bool cond = keepParentEdges, typename = std::enable_if_t<cond>>
  std::vector<int> getReverseEdgePath(int t) {
    assert(distanceLabels[t] != INFTY);
    std::vector<int> path;

    while (parent[t].vertex != INVALID_VERTEX) {
      assert(graph.containsEdge(parent[t].vertex, t));
      assert(graph.edgeHead(parent[t].edge) == t);
      path.push_back(parent[t].edge);
      t = parent[t].vertex;
    }
    return path;
  }

 private:
  // This struct keeps the parent information for a particular vertex v. Depending on the template
  // arguments of the Dijkstra class, we maintain for each vertex only its parent vertex, both its
  // parent vertex and edge, or nothing at all.

  struct ParentVertex {
    int vertex;
  };

  struct FullParentInfo {
    int vertex;
    int edge;
  };

  using ParentInfo = std::conditional_t<keepParentEdges, FullParentInfo, ParentVertex>;

  // Sets the parent vertex of v to u, but only in case we actually keep track of parent vertices.

  template <bool cond = keepParentVertices || keepParentEdges>
  void setParentVertex(const std::enable_if_t<cond, int> v, const int u) {
    assert(v >= 0); assert(v < parent.size());
    parent[v].vertex = u;
  }

  template <bool cond = keepParentVertices || keepParentEdges>
  void setParentVertex(const std::enable_if_t<!cond, int> /*v*/, const int /*u*/) {}

  // Sets the parent edge of v to e, but only in case we actually keep track of parent edges.

  template <bool cond = keepParentEdges>
  void setParentEdge(const std::enable_if_t<cond, int> v, const int e) {
    assert(v >= 0); assert(v < parent.size());
    parent[v].edge = e;
  }

  template <bool cond = keepParentEdges>
  void setParentEdge(const std::enable_if_t<!cond, int> /*v*/, const int /*e*/) {}

  const GraphT& graph;                    // The graph on which we compute shortest paths.
  DistanceLabelContainerT distanceLabels; // The distance labels of the vertices.
  std::vector<ParentInfo> parent;         // The parent information for each vertex.
  QueueT queue;                           // The priority queue of unsettled vertices.
  GetWeightT<GraphT> getWeight;           // A functor returning the edge weight used for routing.
};

// An alias template for a standard Dijkstra search.
template <
    typename GraphT, template <typename> class GetWeightT,
    bool keepParentVertices, bool keepParentEdges>
using StandardDijkstra = Dijkstra<
    GraphT, StampedDistanceLabelContainer, QuadHeap, GetWeightT,
    keepParentVertices, keepParentEdges>;
