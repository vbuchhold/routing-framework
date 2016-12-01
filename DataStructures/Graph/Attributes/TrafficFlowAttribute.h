#pragma once

#include <cassert>

#include "DataStructures/Graph/Attributes/AbstractAttribute.h"

// An attribute associating a traffic flow with each edge of a graph.
class TrafficFlowAttribute : public AbstractAttribute<int> {
 public:
  // A functor that returns the traffic flow of the specified edge in the specified graph. Used for
  // telling algorithms on which attribute of a graph they should work.
  template <typename GraphT>
  struct GetTrafficFlow {
    Type operator()(const GraphT& graph, const int e) const { return graph.trafficFlow(e); }
    Type& operator()(GraphT& graph, const int e) const { return graph.trafficFlow(e); }
  };

  static constexpr const char* NAME = "traffic_flow"; // The attribute's unique name.

  // Returns the traffic flow of edge e.
  Type trafficFlow(const int e) const {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

  // Returns a reference to the traffic flow of edge e.
  Type& trafficFlow(const int e) {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }
};
