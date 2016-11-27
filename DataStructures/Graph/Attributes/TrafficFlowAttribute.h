#pragma once

#include <cassert>

#include "DataStructures/Graph/Attributes/AbstractAttribute.h"

// An attribute associating a traffic flow with each edge of a graph.
class TrafficFlowAttribute : public AbstractAttribute<int> {
 public:
  static constexpr const char* NAME = "traffic_flow"; // The attribute's unique name.

  // Returns the traffic flow of edge e.
  Type trafficFlow(const int e) const {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

  // Sets the traffic flow of edge e to the specified value.
  void setTrafficFlow(const int e, const Type val) {
    assert(e >= 0); assert(e < values.size());
    assert(val >= 0);
    values[e] = val;
  }
};
