#pragma once

#include <cassert>

#include "DataStructures/Graph/Attributes/AbstractAttribute.h"

// An attribute associating a capacity with each edge of a graph.
class CapacityAttribute : public AbstractAttribute<int> {
 public:
  static constexpr const char* NAME = "capacity"; // The attribute's unique name.

  // Returns the capacity in vehicles/h of edge e.
  Type capacity(const int e) const {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

  // Sets the capacity of edge e to the specified value in vehicles/h.
  void setCapacity(const int e, const Type val) {
    assert(e >= 0); assert(e < values.size());
    assert(val >= 0);
    values[e] = val;
  }
};
