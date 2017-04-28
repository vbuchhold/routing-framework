#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include <vectorclass/vectorclass.h>

#include "Algorithms/TrafficAssignment/AllOrNothingAssignment.h"
#include "Algorithms/TrafficAssignment/UnivariateMinimization.h"
#include "DataStructures/Graph/Attributes/TravelCostAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Utilities/OriginDestination.h"
#include "Stats/TrafficAssignment/FrankWolfeAssignmentStats.h"
#include "Tools/Simd/AlignVector.h"
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
        verbose(verbose) {
    stats.totalRunningTime = allOrNothingAssignment.stats.totalRoutingTime;
  }

  // Assigns all OD-flows onto the input graph.
  void run() {
    const AllOrNothingAssignmentStats& substats = allOrNothingAssignment.stats;

    // Initialization.
    Timer timer;
#ifdef TA_NO_SIMD_LINE_SEARCH
    FORALL_EDGES(inputGraph, e)
      inputGraph.travelCost(e) = std::round(objFunction.getEdgeWeight(e, 0));
#else
    FORALL_EDGES_SIMD(inputGraph, e, Vec8f::size()) {
      const Vec8i weight = round_to_int(objFunction.getEdgeWeights(e, 0));
      if (inputGraph.numEdges() - e >= Vec8f::size())
        weight.store(&inputGraph.travelCost(e));
      else
        weight.store_partial(inputGraph.numEdges() - e, &inputGraph.travelCost(e));
    }
#endif
    allOrNothingAssignment.run();
#ifdef TA_NO_SIMD_LINE_SEARCH
    FORALL_EDGES(inputGraph, e) {
      trafficFlows[e] = allOrNothingAssignment.trafficFlowOn(e);
      stats.totalTravelCost += trafficFlows[e] * travelCostFunction(e, trafficFlows[e]);
    }
#else
    Vec8f totalCost = 0;
    FORALL_EDGES_SIMD(inputGraph, e, Vec8f::size()) {
      const Vec8f flow = to_float(Vec8i().load(&allOrNothingAssignment.trafficFlowOn(e)));
      Vec8f cost = flow * travelCostFunction(e, flow);
      if (inputGraph.numEdges() - e >= Vec8f::size()) {
        flow.store(&trafficFlows[e]);
      } else {
        flow.store_partial(inputGraph.numEdges() - e, &trafficFlows[e]);
        cost.cutoff(inputGraph.numEdges() - e);
      }
      totalCost += cost;
    }
    stats.totalTravelCost = horizontal_add(totalCost);
#endif
    stats.lastRunningTime = timer.elapsed();
    stats.lastLineSearchTime = stats.lastRunningTime - substats.lastRoutingTime;
    stats.finishIteration();

    if (verbose) {
      std::cout << "  Line search: " << stats.lastLineSearchTime << "ms";
      std::cout << "  Total: " << stats.lastRunningTime << "ms\n";
      std::cout << "  Total travel cost: " << stats.totalTravelCost << "\n";
      std::cout << std::flush;
    }

    do {
      stats.startIteration();
      Timer timer;

      // Update travel costs.
#ifdef TA_NO_SIMD_LINE_SEARCH
      FORALL_EDGES(inputGraph, e)
        inputGraph.travelCost(e) = std::round(objFunction.getEdgeWeight(e, trafficFlows[e]));
#else
      FORALL_EDGES_SIMD(inputGraph, e, Vec8f::size()) {
        const Vec8f flow = Vec8f().load(&trafficFlows[e]);
        const Vec8i weight = round_to_int(objFunction.getEdgeWeights(e, flow));
        if (inputGraph.numEdges() - e >= Vec8f::size())
          weight.store(&inputGraph.travelCost(e));
        else
          weight.store_partial(inputGraph.numEdges() - e, &inputGraph.travelCost(e));
      }
#endif

      // Direction finding.
      allOrNothingAssignment.run();

      // Line search.
      const float alpha = bisectionMethod([this](const float alpha) {
#ifdef TA_NO_SIMD_LINE_SEARCH
        float sum = 0;
        FORALL_EDGES(inputGraph, e) {
          const float direction = allOrNothingAssignment.trafficFlowOn(e) - trafficFlows[e];
          sum += direction * objFunction.getEdgeWeight(e, trafficFlows[e] + alpha * direction);
        }
        return sum;
#else
        Vec8f sum = 0;
        FORALL_EDGES_SIMD(inputGraph, e, Vec8f::size()) {
          const Vec8f oldFlow = Vec8f().load(&trafficFlows[e]);
          const Vec8f newFlow = to_float(Vec8i().load(&allOrNothingAssignment.trafficFlowOn(e)));
          const Vec8f direction = newFlow - oldFlow;
          Vec8f tmp = direction * objFunction.getEdgeWeights(e, oldFlow + alpha * direction);
          if (inputGraph.numEdges() - e < Vec8f::size())
            tmp.cutoff(inputGraph.numEdges() - e);
          sum += tmp;
        }
        return horizontal_add(sum);
#endif
      }, 0, 1);

      // Move along the descent direction.
#ifdef TA_NO_SIMD_LINE_SEARCH
      FORALL_EDGES(inputGraph, e) {
        const float direction = allOrNothingAssignment.trafficFlowOn(e) - trafficFlows[e];
        trafficFlows[e] = trafficFlows[e] + alpha * direction;
        stats.totalTravelCost += trafficFlows[e] * travelCostFunction(e, trafficFlows[e]);
      }
#else
      Vec8f totalCost = 0;
      FORALL_EDGES_SIMD(inputGraph, e, Vec8f::size()) {
        const Vec8f oldFlow = Vec8f().load(&trafficFlows[e]);
        const Vec8f auxFlow = to_float(Vec8i().load(&allOrNothingAssignment.trafficFlowOn(e)));
        const Vec8f newFlow = oldFlow + alpha * (auxFlow - oldFlow);
        Vec8f cost = newFlow * travelCostFunction(e, newFlow);
        if (inputGraph.numEdges() - e >= Vec8f::size()) {
          newFlow.store(&trafficFlows[e]);
        } else {
          newFlow.store_partial(inputGraph.numEdges() - e, &trafficFlows[e]);
          cost.cutoff(inputGraph.numEdges() - e);
        }
        totalCost += cost;
      }
      stats.totalTravelCost = horizontal_add(totalCost);
#endif
      stats.lastRunningTime = timer.elapsed();
      stats.lastLineSearchTime = stats.lastRunningTime - substats.lastRoutingTime;
      stats.finishIteration();

      if (verbose) {
        std::cout << "  Line search: " << stats.lastLineSearchTime << "ms";
        std::cout << "  Total: " << stats.lastRunningTime << "ms\n";
        std::cout << "  Max change in OD-distances: " << substats.maxChangeInDistances << "\n";
        std::cout << "  Avg change in OD-distances: " << substats.avgChangeInDistances << "\n";
        std::cout << "  Total travel cost: " << stats.totalTravelCost << "\n";
        std::cout << std::flush;
      }
    } while (substats.maxChangeInDistances > 1e-3);

    if (verbose) {
      std::cout << "Total:\n";
      std::cout << "  Checksum: " << substats.totalChecksum;
      std::cout << "  Prepro: " << substats.totalPreprocessingTime << "ms";
      std::cout << "  Custom: " << substats.totalCustomizationTime << "ms";
      std::cout << "  Queries: " << substats.totalQueryTime << "ms";
      std::cout << "  Routing: " << substats.totalRoutingTime << "ms\n";
      std::cout << "  Line search: " << stats.totalLineSearchTime << "ms";
      std::cout << "  Total: " << stats.totalRunningTime << "ms\n";
      std::cout << std::flush;
    }
  }

  // Returns the traffic flow on edge e.
  const float& trafficFlowOn(const int e) const {
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
  AlignVector<float> trafficFlows;       // The traffic flows on the edges.
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
