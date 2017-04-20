#pragma once

#include <cmath>

// The BPR travel cost function, relating the travel time on an edge to the flow on this edge.
template <typename GraphT>
class BprFunction {
 public:
  // Constructs a BPR function.
  BprFunction(const GraphT& graph) : graph(graph) {}

  // Returns the travel time on edge e, given the traffic flow x on e.
  float operator()(const int e, const float x) const {
    return graph.travelTime(e) * (1 + 0.15 * std::pow(x / graph.capacity(e), 4));
  }

  // Returns the derivative of e's travel cost function at x.
  float derivative(const int e, const float x) const {
    const int capacity = graph.capacity(e);
    return graph.travelTime(e) * 0.15 * 4 * std::pow(x / capacity, 4 - 1) / capacity;
  }

 private:
  const GraphT& graph; // The graph on whose edges we operate.
};
