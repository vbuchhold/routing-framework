#pragma once

#include <vectorclass/vectorclass.h>

#include "Tools/Simd/AlignVector.h"

// Represents the system-optimum (SO) objective function. The flow pattern that minimizes the SO
// objective function (while satisfying the flow conservation constraint) minimizes the total
// travel cost. The SO flow pattern is obtained by iterative shortest-path computations using
// appropriate edge weights.
template <typename TravelCostFunctionT>
class SystemOptimum {
 public:
  // Constructs a SO objective function.
  SystemOptimum(TravelCostFunctionT function) : travelCostFunction(function) {}

  // Returns the value of the objective function for the specified edge flows.
  double operator()(const AlignVector<double>& flows) const {
    double sum = 0;
    for (int e = 0; e < flows.size(); ++e)
      sum += flows[e] * travelCostFunction(e, flows[e]);
    return sum;
  }

  // Returns the weight of edge e, given the flow x on e.
  double getEdgeWeight(const int e, const double x) const {
    return travelCostFunction(e, x) + x * travelCostFunction.derivative(e, x);
  }

  // Returns the weights of four consecutive edges starting at e, given the flows x on them.
  Vec4d getEdgeWeights(const int e, const Vec4d& x) const {
    return travelCostFunction(e, x) + x * travelCostFunction.derivative(e, x);
  }

 private:
  TravelCostFunctionT travelCostFunction; // A functor returning the travel cost on an edge.
};
