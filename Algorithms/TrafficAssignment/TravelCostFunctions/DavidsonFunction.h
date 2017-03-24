#pragma once

#include <cassert>
#include <cmath>

// The Davidson travel cost function, relating the travel time on an edge to the flow on this edge.
template <typename GraphT>
class DavidsonFunction {
 public:
  // Constructs a Davidson function.
  DavidsonFunction(const GraphT& graph) : graph(graph) {}

  // Returns the travel time on edge e, given the traffic flow x on e.
  double operator()(const int e, const int x) const {
    assert(x >= 0); assert(x < graph.capacity(e));
    return graph.travelTime(e) * (1 + 0.25 * x / (graph.capacity(e) - x));
  }

  // Returns the derivative of e's travel cost function at x.
  double derivative(const int e, const int x) const {
    assert(x >= 0); assert(x < graph.capacity(e));
    return graph.travelTime(e) * 0.25 * graph.capacity(e) / std::pow(graph.capacity(e) - x, 2);
  }

 private:
  const GraphT& graph; // The graph on whose edges we operate.
};
