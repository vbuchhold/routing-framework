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
  double operator()(const int e, const double x) const {
    const double pt = 0.95 * graph.capacity(e); // The point at which we linearize.
    if (x <= pt)
      return davidson(e, x);
    else
      return davidson(e, pt) + davidson.derivative(e, pt) * (x - pt);
  }

  // Returns the derivative of e's travel cost function at x.
  double derivative(const int e, const double x) const {
    const double pt = 0.95 * graph.capacity(e); // The point at which we linearize.
    if (x <= pt)
      return davidson.derivative(e, x);
    else
      return davidson.derivative(e, pt);
  }

  // Returns the integral of e's travel cost function from 0 to b.
  double integral(const int e, const double b) const {
    const double pt = 0.95 * graph.capacity(e); // The point at which we linearize.
    if (b <= pt)
      return davidson.integral(e, b);
    else
      return davidson.integral(e, pt) + (b - pt) * (operator()(e, b) + operator()(e, pt)) / 2;
  }

  // Returns the travel times on four consecutive edges starting at e, given the flows x on them.
  Vec4d operator()(const int e, const Vec4d& x) const {
    Vec4d pt = 0.95 * to_double(Vec4i().load(&graph.capacity(e)));
    return davidson(e, min(x, pt)) +
        davidson.derivative(e, pt) * (max(x, pt) - pt);
  }

  // Returns the derivative of e's travel cost function at x.
  Vec4d derivative(const int e, const Vec4d& x) const {
    Vec4d pt = 0.95 * to_double(Vec4i().load(&graph.capacity(e)));
    return davidson.derivative(e, min(x, pt));
  }

 private:
  const GraphT& graph;               // The graph on whose edges we operate.
  DavidsonFunction<GraphT> davidson; // The original Davidson function.
};
