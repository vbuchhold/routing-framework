#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
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
  using InputGraph = InputGraphT;

  // Constructs an assignment procedure based on the Frank-Wolfe method.
  FrankWolfeAssignment(InputGraphT& graph, const std::vector<ClusteredOriginDestination>& odPairs,
                       std::ofstream& csv, const bool verbose = true)
      : allOrNothingAssignment(graph, odPairs, verbose),
        inputGraph(graph),
        trafficFlows(graph.numEdges()),
        travelCostFunction(graph),
        objFunction(travelCostFunction),
        csv(csv),
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
    FORALL_EDGES_SIMD(inputGraph, e, Vec4d::size()) {
      const Vec4i weight = round_to_int(objFunction.getEdgeWeights(e, 0));
      if (inputGraph.numEdges() - e >= Vec4d::size())
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
    Vec4d totalCost = 0;
    FORALL_EDGES_SIMD(inputGraph, e, Vec4d::size()) {
      const Vec4d flow = to_double(Vec4i().load(&allOrNothingAssignment.trafficFlowOn(e)));
      Vec4d cost = flow * travelCostFunction(e, flow);
      if (inputGraph.numEdges() - e >= Vec4d::size()) {
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

    if (csv.is_open()) {
      csv << substats.numIterations << "," << substats.lastCustomizationTime << ",";
      csv << substats.lastQueryTime << "," << stats.lastLineSearchTime << ",";
      csv << stats.lastRunningTime << ",nan,nan," << stats.totalTravelCost << ",";
      csv << substats.lastChecksum << std::endl;
    }

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
      FORALL_EDGES_SIMD(inputGraph, e, Vec4d::size()) {
        const Vec4d flow = Vec4d().load(&trafficFlows[e]);
        const Vec4i weight = round_to_int(objFunction.getEdgeWeights(e, flow));
        if (inputGraph.numEdges() - e >= Vec4d::size())
          weight.store(&inputGraph.travelCost(e));
        else
          weight.store_partial(inputGraph.numEdges() - e, &inputGraph.travelCost(e));
      }
#endif

      // Direction finding.
      allOrNothingAssignment.run();

      // Line search.
      const double alpha = bisectionMethod([this](const double alpha) {
#ifdef TA_NO_SIMD_LINE_SEARCH
        double sum = 0;
        FORALL_EDGES(inputGraph, e) {
          const double direction = allOrNothingAssignment.trafficFlowOn(e) - trafficFlows[e];
          sum += direction * objFunction.getEdgeWeight(e, trafficFlows[e] + alpha * direction);
        }
        return sum;
#else
        Vec4d sum = 0;
        FORALL_EDGES_SIMD(inputGraph, e, Vec4d::size()) {
          const Vec4d oldFlow = Vec4d().load(&trafficFlows[e]);
          const Vec4d newFlow = to_double(Vec4i().load(&allOrNothingAssignment.trafficFlowOn(e)));
          const Vec4d direction = newFlow - oldFlow;
          Vec4d tmp = direction * objFunction.getEdgeWeights(e, oldFlow + alpha * direction);
          if (inputGraph.numEdges() - e < Vec4d::size())
            tmp.cutoff(inputGraph.numEdges() - e);
          sum += tmp;
        }
        return horizontal_add(sum);
#endif
      }, 0, 1);

      // Move along the descent direction.
#ifdef TA_NO_SIMD_LINE_SEARCH
      FORALL_EDGES(inputGraph, e) {
        const double direction = allOrNothingAssignment.trafficFlowOn(e) - trafficFlows[e];
        trafficFlows[e] = trafficFlows[e] + alpha * direction;
        stats.totalTravelCost += trafficFlows[e] * travelCostFunction(e, trafficFlows[e]);
      }
#else
      Vec4d totalCost = 0;
      FORALL_EDGES_SIMD(inputGraph, e, Vec4d::size()) {
        const Vec4d oldFlow = Vec4d().load(&trafficFlows[e]);
        const Vec4d auxFlow = to_double(Vec4i().load(&allOrNothingAssignment.trafficFlowOn(e)));
        const Vec4d newFlow = oldFlow + alpha * (auxFlow - oldFlow);
        Vec4d cost = newFlow * travelCostFunction(e, newFlow);
        if (inputGraph.numEdges() - e >= Vec4d::size()) {
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

      if (csv.is_open()) {
        csv << substats.numIterations << "," << substats.lastCustomizationTime << ",";
        csv << substats.lastQueryTime << "," << stats.lastLineSearchTime << ",";
        csv << stats.lastRunningTime << "," << substats.avgChangeInDistances << ",";
        csv << substats.maxChangeInDistances << "," << stats.totalTravelCost << ",";
        csv << substats.lastChecksum << std::endl;
      }

      if (verbose) {
        std::cout << "  Line search: " << stats.lastLineSearchTime << "ms";
        std::cout << "  Total: " << stats.lastRunningTime << "ms\n";
        std::cout << "  Max change in OD-distances: " << substats.maxChangeInDistances << "\n";
        std::cout << "  Avg change in OD-distances: " << substats.avgChangeInDistances << "\n";
        std::cout << "  Total travel cost: " << stats.totalTravelCost << "\n";
        std::cout << std::flush;
      }
    } while (substats.avgChangeInDistances > 1e-2);

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
  const double& trafficFlowOn(const int e) const {
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
  AlignVector<double> trafficFlows;      // The traffic flows on the edges.
  TravelCostFunction travelCostFunction; // A functor returning the travel cost on an edge.
  ObjFunction objFunction;               // The objective function to be minimized (UE or SO).
  std::ofstream& csv;                    // The output CSV file containing statistics.
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
