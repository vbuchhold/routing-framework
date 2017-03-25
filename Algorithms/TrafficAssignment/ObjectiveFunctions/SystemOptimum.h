#pragma once

#include <cmath>

// Represents the system-optimum (SO) objective function. The flow pattern that minimizes the SO
// objective function (while satisfying the flow conservation constraint) minimizes the total
// travel cost. The SO flow pattern is obtained by iterative shortest-path computations using
// appropriate edge weights.
template <typename TravelCostFunctionT>
class SystemOptimum {
 public:
  // Constructs a SO objective function.
  SystemOptimum(TravelCostFunctionT function) : travelCostFunction(function) {}

  // Returns the weight of edge e, given the traffic flow x on e.
  int getEdgeWeight(const int e, const double x) const {
    return std::round(travelCostFunction(e, x) + x * travelCostFunction.derivative(e, x));
  }

 private:
  TravelCostFunctionT travelCostFunction; // A functor returning the travel cost on an edge.
};
