#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "Algorithms/TrafficAssignment/AllOrNothingAssignment.h"
#include "Algorithms/TrafficAssignment/UnivariateMinimization.h"
#include "DataStructures/Graph/Attributes/TravelCostAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Utilities/OriginDestination.h"
#include "Stats/TrafficAssignment/FrankWolfeAssignmentStats.h"
#include "Tools/Timer.h"

// A traffic assignment procedure based on the Frank-Wolfe method (also known as convex combinations
// method). At its heart are iterative shortest-paths computations. The algo can be parameterized to
// compute the user equilibrium or system optimum, and to use different travel cost functions and
// shortest-path algorithms.
template <
    template <typename> class ObjFunctionT, template <typename> class TravelCostFunctionT,
    template <typename, typename> class ShortestPathAlgoT, typename InputGraphT>
class FrankWolfeAssignment {
 public:
  // Constructs an assignment procedure based on the Frank-Wolfe method.
  FrankWolfeAssignment(InputGraphT& graph, const std::vector<OriginDestination>& odPairs,
                       const bool verbose = true)
      : allOrNothingAssignment(graph, odPairs, verbose),
        inputGraph(graph),
        trafficFlows(graph.numEdges()),
        travelCostFunction(graph),
        objFunction(travelCostFunction),
        verbose(verbose) {}

  // Assigns all OD-flows onto the input graph.
  void run() {
    const AllOrNothingAssignmentStats& substats = allOrNothingAssignment.stats;

    // Initialization.
    Timer timer;
    FORALL_EDGES(inputGraph, e)
      inputGraph.travelCost(e) = objFunction.getEdgeWeight(e, 0);
    allOrNothingAssignment.run();
    FORALL_EDGES(inputGraph, e) {
      trafficFlows[e] = allOrNothingAssignment.trafficFlowOn(e);
      stats.totalTravelCost += trafficFlows[e] * travelCostFunction(e, trafficFlows[e]);
    }
    stats.lastRunningTime = timer.elapsed();
    stats.finishIteration();

    if (verbose) {
      std::cout << "  Running time: " << stats.lastRunningTime << "ms\n";
      std::cout << "  Total travel cost: " << stats.totalTravelCost << "\n";
      std::cout << std::flush;
    }

    do {
      stats.startIteration();
      Timer timer;

      // Update travel costs.
      FORALL_EDGES(inputGraph, e)
        inputGraph.travelCost(e) = objFunction.getEdgeWeight(e, trafficFlows[e]);

      // Direction finding.
      allOrNothingAssignment.run();

      // Line search.
      const float alpha = bisectionMethod([this](const float alpha) {
        float sum = 0;
        FORALL_EDGES(inputGraph, e) {
          const float direction = allOrNothingAssignment.trafficFlowOn(e) - trafficFlows[e];
          sum += direction * objFunction.getEdgeWeight(e, trafficFlows[e] + alpha * direction);
        }
        return sum;
      }, 0, 1);

      // Move along the descent direction.
      float totalFlow = 0;
      FORALL_EDGES(inputGraph, e) {
        const float direction = allOrNothingAssignment.trafficFlowOn(e) - trafficFlows[e];
        stats.changeInEdgeFlows += std::abs(alpha * direction);
        totalFlow += trafficFlows[e];
        trafficFlows[e] = trafficFlows[e] + alpha * direction;
        stats.totalTravelCost += trafficFlows[e] * travelCostFunction(e, trafficFlows[e]);
      }
      stats.lastRunningTime = timer.elapsed();
      stats.changeInEdgeFlows /= totalFlow;
      stats.finishIteration();

      if (verbose) {
        std::cout << "  Running time: " << stats.lastRunningTime << "ms\n";
        std::cout << "  Max change in OD-distances: " << substats.maxChangeInDistances << "\n";
        std::cout << "  Avg change in OD-distances: " << substats.avgChangeInDistances << "\n";
        std::cout << "  Change in edge flows: " << stats.changeInEdgeFlows << "\n";
        std::cout << "  Total travel cost: " << stats.totalTravelCost << "\n";
        std::cout << std::flush;
      }
    } while (substats.maxChangeInDistances > 1e-3);

    if (verbose) {
      std::cout << "Total:\n";
      std::cout << "  Checksum: " << substats.totalChecksum;
      std::cout << "  Prepro: " << substats.totalPreprocessingTime << "ms";
      std::cout << "  Custom: " << substats.totalCustomizationTime << "ms";
      std::cout << "  Queries: " << substats.totalQueryTime << "ms\n";
      std::cout << "  Running time: " << stats.totalRunningTime << "ms\n";
      std::cout << std::flush;
    }
  }

  // Returns the traffic flow on edge e.
  float trafficFlowOn(const int e) const {
    assert(e >= 0); assert(e < inputGraph.numEdges());
    return trafficFlows[e];
  }

  FrankWolfeAssignmentStats stats; // Statistics about the execution.

 private:
  using AllOrNothing = AllOrNothingAssignment<ShortestPathAlgoT<InputGraphT, TravelCostAttribute>>;
  using TravelCostFunction = TravelCostFunctionT<InputGraphT>;
  using ObjFunction = ObjFunctionT<TravelCostFunction>;

  AllOrNothing allOrNothingAssignment;   // The all-or-nothing assignment algo used as a subroutine.
  InputGraphT& inputGraph;               // The input graph.
  std::vector<float> trafficFlows;       // The traffic flows on the edges.
  TravelCostFunction travelCostFunction; // A functor returning the travel cost on an edge.
  ObjFunction objFunction;               // The objective function to be minimized (UE or SO).
  const bool verbose;                    // Should informative messages be displayed?
};

// An alias template for a user-equilibrium (UE) traffic assignment.
template <
    template <typename> class TravelCostFunctionT,
    template <typename, typename> class ShortestPathAlgoT, typename InputGraphT>
using UEAssignment =
    FrankWolfeAssignment<UserEquilibrium, TravelCostFunctionT, ShortestPathAlgoT, InputGraphT>;

// An alias template for a system-optimum (SO) traffic assignment.
template <
    template <typename> class TravelCostFunctionT,
    template <typename, typename> class ShortestPathAlgoT, typename InputGraphT>
using SOAssignment =
    FrankWolfeAssignment<SystemOptimum, TravelCostFunctionT, ShortestPathAlgoT, InputGraphT>;
