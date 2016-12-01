#pragma once

#include <cassert>

#include "DataStructures/Graph/Attributes/AbstractAttribute.h"

// An attribute associating a capacity with each edge of a graph.
class CapacityAttribute : public AbstractAttribute<int> {
 public:
  // A functor that returns the capacity of the specified edge in the specified graph. Used for
  // telling algorithms on which attribute of a graph they should work.
  template <typename GraphT>
  struct GetCapacity {
    Type operator()(const GraphT& graph, const int e) const { return graph.capacity(e); }
    Type& operator()(GraphT& graph, const int e) const { return graph.capacity(e); }
  };

  static constexpr const char* NAME = "capacity"; // The attribute's unique name.

  // Returns the capacity in vehicles/h of edge e.
  Type capacity(const int e) const {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

  // Returns a reference to the capacity in vehicles/h of edge e.
  Type& capacity(const int e) {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }
};
