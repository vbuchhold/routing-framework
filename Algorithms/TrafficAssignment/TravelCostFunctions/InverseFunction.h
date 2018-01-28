#pragma once

#include <cassert>
#include <cmath>

#include <vectorclass/vectorclass.h>

// An inverse travel cost function, which decreases (rather than increases) as the flow increases.
// It should model the costs for operating public transit. The more people use public transit, the
// lower the costs per person are. However, a traffic assignment minimizing solely the operational
// cost may result in long detours for some passengers. Therefore, the (static) travel time also
// contributes to the travel cost.
template <typename GraphT>
class InverseFunction {
 public:
  // Constructs an inverse travel cost function.
  InverseFunction(const GraphT& graph) : graph(graph) {}

  // Returns the travel cost on edge e, given the flow x on e.
  double operator()(const int e, const double x) const {
    assert(x >= 0);
    return 1.0 * graph.length(e) * (155 / (x + 1) + 0.85) + 0.0 * graph.travelTime(e);
  }

  // Returns the derivative of e's travel cost function at x.
  double derivative(const int e, const double x) const {
    assert(x >= 0);
    const double tmp = x + 1;
    return 1.0 * graph.length(e) * -155 / (tmp * tmp);
  }

  // Returns the antiderivative of e's travel cost function at x.
  double antiderivative(const int e, const double x) const {
    assert(x >= 0);
    const int length = graph.length(e);
    const int time = graph.travelTime(e);
    return 1.0 * length * (155 * std::log(x + 1) + 0.85 * x) + 0.0 * time * x;
  }

  // Returns the integral of e's travel cost function from 0 to b.
  double integral(const int e, const double b) const {
    return antiderivative(e, b) - antiderivative(e, 0);
  }

  // Returns the travel costs on four consecutive edges starting at e, given the flows x on them.
  Vec4d operator()(const int e, const Vec4d& x) const {
    Vec4d length = to_double(Vec4i().load(&graph.length(e)));
    Vec4d time = to_double(Vec4i().load(&graph.travelTime(e)));
    return 1.0 * length * (155 / (x + 1) + 0.85) + 0.0 * time;
  }

  // Returns the derivative of e's travel cost function at x.
  Vec4d derivative(const int e, const Vec4d& x) const {
    Vec4d length = to_double(Vec4i().load(&graph.length(e)));
    return 1.0 * length * -155 / pow_const(x + 1, 2);
  }

 private:
  const GraphT& graph; // The graph on whose edges we operate.
};
