#pragma once

#include <cassert>

#include "DataStructures/Graph/Attributes/AbstractAttribute.h"

// An attribute associating a capacity with each edge of a graph.
class CapacityAttribute : public AbstractAttribute<int> {
 public:
  // Returns the capacity in vehicles/h of edge e.
  const Type& capacity(const int e) const {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

  // Returns a reference to the capacity in vehicles/h of edge e.
  Type& capacity(const int e) {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

 protected:
  static constexpr const char* NAME = "capacity"; // The attribute's unique name.
};
