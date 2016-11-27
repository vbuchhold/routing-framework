#pragma once

#include <cassert>

#include "DataStructures/Graph/Attributes/AbstractAttribute.h"

// An attribute associating a routing cost with each edge of a graph.
class RoutingCostAttribute : public AbstractAttribute<int> {
 public:
  static constexpr Type DEFAULT_VALUE = INFTY;          // The attribute's default value.
  static constexpr const char* NAME   = "routing_cost"; // The attribute's unique name.

  // Returns the routing cost of edge e.
  Type routingCost(const int e) const {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

  // Sets the routing cost of edge e to the specified value.
  void setRoutingCost(const int e, const Type val) {
    assert(e >= 0); assert(e < values.size());
    values[e] = val;
  }
};
