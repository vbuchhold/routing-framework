#pragma once

#include <cassert>

#include "DataStructures/Graph/Attributes/AbstractAttribute.h"

// An attribute associating an ID with each edge of a graph.
class EdgeIdAttribute : public AbstractAttribute<int> {
 public:
  // A functor that returns the ID of the specified edge in the specified graph. Used for telling
  // algorithms on which attribute of a graph they should work.
  template <typename GraphT>
  struct GetEdgeId {
    Type operator()(const GraphT& graph, const int e) const { return graph.edgeId(e); }
    Type& operator()(GraphT& graph, const int e) const { return graph.edgeId(e); }
  };

  static constexpr Type DEFAULT_VALUE = -1;        // The attribute's default value.
  static constexpr const char* NAME   = "edge_id"; // The attribute's unique name.

  // Returns the ID of edge e.
  Type edgeId(const int e) const {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

  // Returns a reference to the ID of edge e.
  Type& edgeId(const int e) {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }
};
