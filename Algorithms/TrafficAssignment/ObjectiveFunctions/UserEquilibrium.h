#pragma once

#include <vector>

// Represents the user-equilibrium (UE) objective function. The flow pattern that minimizes the UE
// objective function (while satisfying the flow conservation constraint) is such that all drivers
// minimize their own travel cost. The UE flow pattern is obtained by iterative shortest-path
// computations using appropriate edge weights.
template <typename TravelCostFunctionT>
class UserEquilibrium {
 public:
  // Constructs an UE objective function.
  UserEquilibrium(TravelCostFunctionT function) : travelCostFunction(function) {}

  // Returns the value of the objective function for the specified edge flows.
  double operator()(const std::vector<double>& flows) const {
    double sum = 0;
    for (int e = 0; e < flows.size(); ++e)
      sum += travelCostFunction.integral(e, flows[e]);
    return sum;
  }

  // Returns the partial derivative with respect to the e-th variable x_e at x_e = x.
  double derivative(const int e, const double x) const {
    return travelCostFunction(e, x);
  }

  // Returns the second order partial derivative with respect to the e-th variable x_e at x_e = x.
  double secondDerivative(const int e, const double x) const {
    return travelCostFunction.derivative(e, x);
  }

 private:
  TravelCostFunctionT travelCostFunction; // A functor returning the travel cost on an edge.
};
