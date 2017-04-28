#pragma once

#include <vectorclass/vectorclass.h>

#include "Algorithms/TrafficAssignment/TravelCostFunctions/DavidsonFunction.h"

// The modified Davidson travel cost function. It relates the travel time on an edge to the flow on
// this edge.
template <typename GraphT>
class ModifiedDavidsonFunction {
 public:
  // Constructs a modified Davidson function.
  ModifiedDavidsonFunction(const GraphT& graph) : graph(graph), davidson(graph) {}

  // Returns the travel time on edge e, given the flow x on e.
  float operator()(const int e, const float x) const {
    const float operatingPt = 0.95f * graph.capacity(e); // The point at which we linearize.
    if (x <= operatingPt)
      return davidson(e, x);
    else
      return davidson(e, operatingPt) + davidson.derivative(e, operatingPt) * (x - operatingPt);
  }

  // Returns the derivative of e's travel cost function at x.
  float derivative(const int e, const float x) const {
    const float operatingPt = 0.95f * graph.capacity(e); // The point at which we linearize.
    if (x <= operatingPt)
      return davidson.derivative(e, x);
    else
      return davidson.derivative(e, operatingPt);
  }

  // Returns the travel times on eight consecutive edges starting at e, given the flows x on them.
  Vec8f operator()(const int e, const Vec8f& x) const {
    Vec8f operatingPt = 0.95f * to_float(Vec8i().load(&graph.capacity(e)));
    return davidson(e, min(x, operatingPt)) +
        davidson.derivative(e, operatingPt) * (max(x, operatingPt) - operatingPt);
  }

  // Returns the derivative of e's travel cost function at x.
  Vec8f derivative(const int e, const Vec8f& x) const {
    Vec8f operatingPt = 0.95f * to_float(Vec8i().load(&graph.capacity(e)));
    return davidson.derivative(e, min(x, operatingPt));
  }

 private:
  const GraphT& graph;               // The graph on whose edges we operate.
  DavidsonFunction<GraphT> davidson; // The original Davidson function.
};
