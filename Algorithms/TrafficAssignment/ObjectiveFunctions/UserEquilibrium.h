#pragma once

#include <vectorclass/vectorclass.h>

// Represents the user-equilibrium (UE) objective function. The flow pattern that minimizes the UE
// objective function (while satisfying the flow conservation constraint) is such that all drivers
// minimize their own travel cost. The UE flow pattern is obtained by iterative shortest-path
// computations using appropriate edge weights.
template <typename TravelCostFunctionT>
class UserEquilibrium {
 public:
  // Constructs an UE objective function.
  UserEquilibrium(TravelCostFunctionT function) : travelCostFunction(function) {}

  // Returns the weight of edge e, given the flow x on e.
  float getEdgeWeight(const int e, const float x) const {
    return travelCostFunction(e, x);
  }

  // Returns the weights of eight consecutive edges starting at e, given the flows x on them.
  Vec8f getEdgeWeights(const int e, const Vec8f& x) const {
    return travelCostFunction(e, x);
  }

 private:
  TravelCostFunctionT travelCostFunction; // A functor returning the travel cost on an edge.
};
