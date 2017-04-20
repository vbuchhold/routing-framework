#pragma once

#include <cassert>

#include <vectorclass/vectorclass.h>

// The Davidson travel cost function, relating the travel time on an edge to the flow on this edge.
template <typename GraphT>
class DavidsonFunction {
 public:
  // Constructs a Davidson function.
  DavidsonFunction(const GraphT& graph) : graph(graph) {}

  // Returns the travel time on edge e, given the flow x on e.
  float operator()(const int e, const float x) const {
    assert(x >= 0); assert(x < graph.capacity(e));
    return graph.travelTime(e) * (1 + 0.25f * x / (graph.capacity(e) - x));
  }

  // Returns the derivative of e's travel cost function at x.
  float derivative(const int e, const float x) const {
    assert(x >= 0); assert(x < graph.capacity(e));
    const float tmp = graph.capacity(e) - x;
    return graph.travelTime(e) * 0.25f * graph.capacity(e) / (tmp * tmp);
  }

  // Returns the travel times on eight consecutive edges starting at e, given the flows x on them.
  Vec8f operator()(const int e, const Vec8f& x) const {
    Vec8f time = to_float(Vec8i().load(&graph.travelTime(e)));
    Vec8f capacity = to_float(Vec8i().load(&graph.capacity(e)));
    return time * (1 + 0.25f * x / (capacity - x));
  }

  // Returns the derivative of e's travel cost function at x.
  Vec8f derivative(const int e, const Vec8f& x) const {
    Vec8f time = to_float(Vec8i().load(&graph.travelTime(e)));
    Vec8f capacity = to_float(Vec8i().load(&graph.capacity(e)));
    return time * 0.25f * capacity / pow_const(capacity - x, 2);
  }

 private:
  const GraphT& graph; // The graph on whose edges we operate.
};
