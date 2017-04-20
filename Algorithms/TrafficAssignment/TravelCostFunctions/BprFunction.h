#pragma once

#include <vectorclass/vectorclass.h>

// The BPR travel cost function, relating the travel time on an edge to the flow on this edge.
template <typename GraphT>
class BprFunction {
 public:
  // Constructs a BPR function.
  BprFunction(const GraphT& graph) : graph(graph) {}

  // Returns the travel time on edge e, given the flow x on e.
  float operator()(const int e, const float x) const {
    const float tmp = x / graph.capacity(e);
    return graph.travelTime(e) * (1 + 0.15f * tmp * tmp * tmp * tmp);
  }

  // Returns the derivative of e's travel cost function at x.
  float derivative(const int e, const float x) const {
    const float tmp = x / graph.capacity(e);
    return graph.travelTime(e) * 0.15f * 4 * tmp * tmp * tmp / graph.capacity(e);
  }

  // Returns the travel times on eight consecutive edges starting at e, given the flows x on them.
  Vec8f operator()(const int e, const Vec8f& x) const {
    Vec8f time = to_float(Vec8i().load(&graph.travelTime(e)));
    Vec8f capacity = to_float(Vec8i().load(&graph.capacity(e)));
    return time * (1 + 0.15f * pow_const(x / capacity, 4));
  }

  // Returns the derivative of e's travel cost function at x.
  Vec8f derivative(const int e, const Vec8f& x) const {
    Vec8f time = to_float(Vec8i().load(&graph.travelTime(e)));
    Vec8f capacity = to_float(Vec8i().load(&graph.capacity(e)));
    return time * 0.15f * 4 * pow_const(x / capacity, 4 - 1) / capacity;
  }

 private:
  const GraphT& graph; // The graph on whose edges we operate.
};
