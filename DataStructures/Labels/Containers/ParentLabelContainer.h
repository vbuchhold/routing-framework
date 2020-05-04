#pragma once

#include <cassert>
#include <cstdint>
#include <type_traits>
#include <vector>

#include "Tools/Simd/AlignedVector.h"
#include "Tools/Constants.h"

// A container that maintains parent information during a shortest-path search. Depending on the
// used label set, the parent information consists of parent vertices and/or edges.
template <typename GraphT, typename LabelSetT>
class ParentLabelContainer {
 public:
  using LabelMask = typename LabelSetT::LabelMask;     // Marks a subset of components in a label.
  using ParentLabel = typename LabelSetT::ParentLabel; // The parent information for a vertex.

  // Constructs a container that maintains parent information.
  explicit ParentLabelContainer(const GraphT& graph)
      : graph(graph), parent(LabelSetT::KEEP_PARENT_VERTICES ? graph.numVertices() : 0) {}

  // Sets the parent vertex of v to u on all shortest paths specified by mask.

  template <bool COND = LabelSetT::KEEP_PARENT_VERTICES>
  std::enable_if_t<COND> setVertex(const int v, const int u, const LabelMask& mask) {
    assert(v >= 0); assert(v < parent.size());
    parent[v].setVertex(u, mask);
  }

  template <bool COND = LabelSetT::KEEP_PARENT_VERTICES>
  std::enable_if_t<!COND> setVertex(const int /*v*/, const int /*u*/, const LabelMask& /*mask*/) {}

  // Sets the parent edge of v to e on all shortest paths specified by mask.

  template <bool COND = LabelSetT::KEEP_PARENT_EDGES>
  std::enable_if_t<COND> setEdge(const int v, const int e, const LabelMask& mask) {
    assert(v >= 0); assert(v < parent.size());
    parent[v].setEdge(e, mask);
  }

  template <bool COND = LabelSetT::KEEP_PARENT_EDGES>
  std::enable_if_t<!COND> setEdge(const int /*v*/, const int /*e*/, const LabelMask& /*mask*/) {}

  // Returns the parent vertex of v on the shortest path from the i-th source to v.
  int getVertex(const int v, const int i = 0) const {
    static_assert(LabelSetT::KEEP_PARENT_VERTICES, "We currently do not keep parent vertices.");
    assert(v >= 0); assert(v < parent.size());
    return parent[v].vertex(i);
  }

  // Returns the parent edge of v on the shortest path from the i-th source to v.
  int getEdge(const int v, const int i = 0) const {
    static_assert(LabelSetT::KEEP_PARENT_EDGES, "We currently do not keep parent edges.");
    assert(v >= 0); assert(v < parent.size());
    return parent[v].edge(i);
  }

  // Returns the vertices on the shortest path from the i-th source to t in reverse order.
  const std::vector<int32_t>& getReversePath(int t, const int i = 0) {
    static_assert(LabelSetT::KEEP_PARENT_VERTICES, "We currently do not keep parent vertices.");
    assert(t >= 0); assert(t < parent.size());
    lastPath.clear();
    while (parent[t].vertex(i) != t) {
      assert(graph.containsEdge(parent[t].vertex(i), t));
      lastPath.push_back(t);
      t = parent[t].vertex(i);
    }
    lastPath.push_back(t);
    return lastPath;
  }

  // Returns the edges on the shortest path from the i-th source to t in reverse order.
  const std::vector<int32_t>& getReverseEdgePath(int t, const int i = 0) {
    static_assert(LabelSetT::KEEP_PARENT_EDGES, "We currently do not keep parent edges.");
    assert(t >= 0); assert(t < parent.size());
    lastPath.clear();
    while (parent[t].vertex(i) != t) {
      assert(graph.containsEdge(parent[t].vertex(i), t));
      assert(graph.edgeHead(parent[t].edge(i)) == t);
      lastPath.push_back(parent[t].edge(i));
      t = parent[t].vertex(i);
    }
    return lastPath;
  }

 private:
  const GraphT& graph;               // The graph on which we compute shortest paths.
  AlignedVector<ParentLabel> parent; // The parent information for each vertex.
  std::vector<int32_t> lastPath;     // The (vertex or edge) path that was retrieved last.
};
