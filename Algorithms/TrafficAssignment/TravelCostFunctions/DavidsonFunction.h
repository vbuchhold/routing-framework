#pragma once

#include <cassert>
#include <cmath>

#include <vectorclass/vectorclass.h>

// The Davidson travel cost function, relating the travel time on an edge to the flow on this edge.
template <typename GraphT>
class DavidsonFunction {
 public:
  // Constructs a Davidson function.
  DavidsonFunction(const GraphT& graph) : graph(graph) {}

  // Returns the travel time on edge e, given the flow x on e.
  double operator()(const int e, const double x) const {
    assert(x >= 0); assert(x < graph.capacity(e));
    return graph.travelTime(e) * (1 + 0.01 * x / (graph.capacity(e) - x));
  }

  // Returns the derivative of e's travel cost function at x.
  double derivative(const int e, const double x) const {
    assert(x >= 0); assert(x < graph.capacity(e));
    const double tmp = graph.capacity(e) - x;
    return graph.travelTime(e) * 0.01 * graph.capacity(e) / (tmp * tmp);
  }

  // Returns the antiderivative of e's travel cost function at x.
  double antiderivative(const int e, const double x) const {
    assert(x >= 0); assert(x < graph.capacity(e));
    const int time = graph.travelTime(e);
    const int capacity = graph.capacity(e);
    return time * (x - 0.01 * (capacity * std::log(capacity - x) + x));
  }

  // Returns the integral of e's travel cost function from 0 to b.
  double integral(const int e, const double b) const {
    return antiderivative(e, b) - antiderivative(e, 0);
  }

  // Returns the travel times on four consecutive edges starting at e, given the flows x on them.
  Vec4d operator()(const int e, const Vec4d& x) const {
    Vec4d time = to_double(Vec4i().load(&graph.travelTime(e)));
    Vec4d capacity = to_double(Vec4i().load(&graph.capacity(e)));
    return time * (1 + 0.01 * x / (capacity - x));
  }

  // Returns the derivative of e's travel cost function at x.
  Vec4d derivative(const int e, const Vec4d& x) const {
    Vec4d time = to_double(Vec4i().load(&graph.travelTime(e)));
    Vec4d capacity = to_double(Vec4i().load(&graph.capacity(e)));
    return time * 0.01 * capacity / pow_const(capacity - x, 2);
  }

 private:
  const GraphT& graph; // The graph on whose edges we operate.
};
