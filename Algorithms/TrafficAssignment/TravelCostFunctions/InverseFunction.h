#pragma once

#include <cassert>

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
  float operator()(const int e, const float x) const {
    assert(x >= 0);
    return 1.0f * graph.length(e) * (155 / (x + 1) + 0.85f) + 0.0f * graph.travelTime(e);
  }

  // Returns the derivative of e's travel cost function at x.
  float derivative(const int e, const float x) const {
    assert(x >= 0);
    const float tmp = x + 1;
    return 1.0f * graph.length(e) * -155 / (tmp * tmp);
  }

  // Returns the travel costs on eight consecutive edges starting at e, given the flows x on them.
  Vec8f operator()(const int e, const Vec8f& x) const {
    Vec8f length = to_float(Vec8i().load(&graph.length(e)));
    Vec8f time = to_float(Vec8i().load(&graph.travelTime(e)));
    return 1.0f * length * (155 / (x + 1) + 0.85f) + 0.0f * time;
  }

  // Returns the derivative of e's travel cost function at x.
  Vec8f derivative(const int e, const Vec8f& x) const {
    Vec8f length = to_float(Vec8i().load(&graph.length(e)));
    return 1.0f * length * -155 / pow_const(x + 1, 2);
  }

 private:
  const GraphT& graph; // The graph on whose edges we operate.
};
