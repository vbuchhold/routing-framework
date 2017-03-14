#pragma once

#include <cassert>

#include "DataStructures/Graph/Attributes/AbstractAttribute.h"
#include "Tools/Constants.h"

// An attribute associating a routing cost with each edge of a graph.
class RoutingCostAttribute : public AbstractAttribute<int> {
 public:
  static constexpr Type DEFAULT_VALUE = INFTY;          // The attribute's default value.
  static constexpr const char* NAME   = "routing_cost"; // The attribute's unique name.

  // Returns the routing cost of edge e.
  const Type& routingCost(const int e) const {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

  // Returns a reference to the routing cost of edge e.
  Type& routingCost(const int e) {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }
};
