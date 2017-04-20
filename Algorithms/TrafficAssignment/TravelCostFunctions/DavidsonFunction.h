#pragma once

#include <cassert>

// The Davidson travel cost function, relating the travel time on an edge to the flow on this edge.
template <typename GraphT>
class DavidsonFunction {
 public:
  // Constructs a Davidson function.
  DavidsonFunction(const GraphT& graph) : graph(graph) {}

  // Returns the travel time on edge e, given the traffic flow x on e.
  float operator()(const int e, const float x) const {
    assert(x >= 0); assert(x < graph.capacity(e));
    return graph.travelTime(e) * (1 + 0.25 * x / (graph.capacity(e) - x));
  }

  // Returns the derivative of e's travel cost function at x.
  float derivative(const int e, const float x) const {
    assert(x >= 0); assert(x < graph.capacity(e));
    const float tmp = graph.capacity(e) - x;
    return graph.travelTime(e) * 0.25 * graph.capacity(e) / (tmp * tmp);
  }

 private:
  const GraphT& graph; // The graph on whose edges we operate.
};
