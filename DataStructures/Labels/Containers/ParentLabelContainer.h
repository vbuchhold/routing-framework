#pragma once

#include <cassert>
#include <type_traits>
#include <vector>

#include "Tools/Simd/AlignedVector.h"

// A container maintaining parent information during a shortest-path search. Depending on the used
// label set, the parent information consists of parent vertices and/or edges.
template <typename GraphT, typename LabelSetT>
class ParentLabelContainer {
 public:
  using LabelMask = typename LabelSetT::LabelMask;     // Marks a subset of components in a label.
  using ParentLabel = typename LabelSetT::ParentLabel; // The parent information for a vertex.

  // Constructs a container maintaining parent information.
  ParentLabelContainer(const GraphT& graph) : graph(graph) {
    resize();
  }

  // Ensures that the container can hold as many parent labels as the graph has vertices.
  void resize() {
    parent.resize(LabelSetT::KEEP_PARENT_VERTICES ? graph.numVertices() : 0);
  }

  // Sets the parent vertex of v to u on all shortest paths specified by mask.

  template <bool cond = LabelSetT::KEEP_PARENT_VERTICES>
  std::enable_if_t<cond> setVertex(const int v, const int u, const LabelMask& mask) {
    assert(v >= 0); assert(v < parent.size());
    parent[v].setVertex(u, mask);
  }

  template <bool cond = LabelSetT::KEEP_PARENT_VERTICES>
  std::enable_if_t<!cond> setVertex(const int /*v*/, const int /*u*/, const LabelMask& /*mask*/) {}

  // Sets the parent edge of v to e on all shortest paths specified by mask.

  template <bool cond = LabelSetT::KEEP_PARENT_EDGES>
  std::enable_if_t<cond> setEdge(const int v, const int e, const LabelMask& mask) {
    assert(v >= 0); assert(v < parent.size());
    parent[v].setEdge(e, mask);
  }

  template <bool cond = LabelSetT::KEEP_PARENT_EDGES>
  std::enable_if_t<!cond> setEdge(const int /*v*/, const int /*e*/, const LabelMask& /*mask*/) {}

  // Returns the vertices on the shortest path from the i-th source to t in reverse order.
  std::vector<int> getReversePath(int t, const int i = 0) const {
    static_assert(LabelSetT::KEEP_PARENT_VERTICES, "We currently do not keep parent vertices.");
    std::vector<int> path;
    while (parent[t].vertex(i) != INVALID_VERTEX) {
      assert(graph.containsEdge(parent[t].vertex(i), t));
      path.push_back(t);
      t = parent[t].vertex(i);
    }
    path.push_back(t);
    return path;
  }

  // Returns the edges on the shortest path from the i-th source to t in reverse order.
  std::vector<int> getReverseEdgePath(int t, const int i = 0) const {
    static_assert(LabelSetT::KEEP_PARENT_EDGES, "We currently do not keep parent edges.");
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
  const GraphT& graph;               // The graph on which we compute shortest paths.
  AlignedVector<ParentLabel> parent; // The parent information for each vertex.
};
