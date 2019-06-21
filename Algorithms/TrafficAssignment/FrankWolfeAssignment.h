#pragma once

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

#include "Algorithms/TrafficAssignment/ObjectiveFunctions/SystemOptimum.h"
#include "Algorithms/TrafficAssignment/ObjectiveFunctions/UserEquilibrium.h"
#include "Algorithms/TrafficAssignment/AllOrNothingAssignment.h"
#include "Algorithms/TrafficAssignment/UnivariateMinimization.h"
#include "DataStructures/Graph/Attributes/TraversalCostAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Utilities/OriginDestination.h"
#include "Stats/TrafficAssignment/FrankWolfeAssignmentStats.h"
#include "Tools/Math.h"
#include "Tools/Timer.h"

// A traffic assignment procedure based on the Frank-Wolfe method (also known as convex combinations
// method). At its heart are iterative shortest-paths computations. The algo can be parameterized to
// compute the user equilibrium or system optimum, and to use different traversal cost functions and
// shortest-path algorithms.
template <
    template <typename> class ObjFunctionT, template <typename> class TraversalCostFunctionT,
    template <typename, typename> class ShortestPathAlgoT, typename GraphT>
class FrankWolfeAssignment {
 public:
  using Graph = GraphT;

  // Constructs an assignment procedure based on the Frank-Wolfe method.
  FrankWolfeAssignment(Graph& graph, const std::vector<ClusteredOriginDestination>& odPairs,
                       const bool verbose = true)
      : aonAssignment(graph, odPairs, verbose),
        graph(graph),
        trafficFlows(graph.numEdges()),
        pointOfSight(graph.numEdges()),
        traversalCostFunction(graph),
        objFunction(traversalCostFunction),
        verbose(verbose) {
    assert(graph.isDefrag());
    stats.totalRunningTime = aonAssignment.stats.totalRoutingTime;
  }

  // Assigns all OD flows onto the graph.
  void run(
      std::ofstream& flowFile, std::ofstream& distFile, std::ofstream& statFile,
      const int numIterations = 0, const bool outputIntermediates = false) {
    assert(numIterations >= 0);
    Timer timer;
    auto prevSkipInterval = 1u;
    determineInitialSolution(prevSkipInterval);
    stats.lastRunningTime = timer.elapsed();
    stats.lastLineSearchTime = stats.lastRunningTime - aonAssignment.stats.lastRoutingTime;
    stats.finishIteration();

    if (flowFile.is_open())
      FORALL_EDGES(graph, e) {
        const auto vol = trafficFlows[e];
        const auto sat = vol / graph.capacity(e);
        flowFile << aonAssignment.stats.numIterations << ',' << vol << ',' << sat << '\n';
      }

    if (distFile.is_open())
      for (const auto dist : aonAssignment.stats.lastDistances)
        distFile << aonAssignment.stats.numIterations << ',' << dist << '\n';

    if (statFile.is_open()) {
      statFile << aonAssignment.stats.numIterations << ",";
      statFile << aonAssignment.stats.lastCustomizationTime << ",";
      statFile << aonAssignment.stats.lastQueryTime << ",";
      statFile << stats.lastLineSearchTime << "," << stats.lastRunningTime << ",nan,nan,";
      statFile << aonAssignment.stats.lastChecksum << std::endl;
    }

    if (verbose) {
      std::cout << "  Line search: " << stats.lastLineSearchTime << "ms";
      std::cout << "  Total: " << stats.lastRunningTime << "ms\n";
      std::cout << std::flush;
    }

    while ((numIterations != 0 || stats.prevRelGap > 1e-4 || prevSkipInterval > 1) &&
           (numIterations == 0 || aonAssignment.stats.numIterations < numIterations)) {
      stats.startIteration();
      Timer timer;
      const unsigned int skip = std::min(std::max(stats.prevRelGap / 1e-4, 1.0), double{-1u});
      const auto skipInterval = std::min(roundDownToPowerOfTwo(skip), prevSkipInterval);
      updateTraversalCosts();
      findDescentDirection(skipInterval);
      const auto tau = findMoveSize();
      moveAlongDescentDirection(tau);
      const auto prevMinPathCost = aonAssignment.stats.prevMinPathCost * prevSkipInterval;
      assert(prevMinPathCost <= stats.prevTotalPathCost);
      prevSkipInterval = skipInterval;
      stats.lastRunningTime = timer.elapsed();
      stats.lastLineSearchTime = stats.lastRunningTime - aonAssignment.stats.lastRoutingTime;
      stats.prevRelGap = 1 - prevMinPathCost / stats.prevTotalPathCost;
      stats.finishIteration();

      if (flowFile.is_open() && outputIntermediates)
        FORALL_EDGES(graph, e) {
          const auto vol = trafficFlows[e];
          const auto sat = vol / graph.capacity(e);
          flowFile << aonAssignment.stats.numIterations << ',' << vol << ',' << sat << '\n';
        }

      if (distFile.is_open() && outputIntermediates)
        for (const auto dist : aonAssignment.stats.lastDistances)
          distFile << aonAssignment.stats.numIterations << ',' << dist << '\n';

      if (statFile.is_open()) {
        statFile << aonAssignment.stats.numIterations << ",";
        statFile << aonAssignment.stats.lastCustomizationTime << ",";
        statFile << aonAssignment.stats.lastQueryTime << ",";
        statFile << stats.lastLineSearchTime << "," << stats.lastRunningTime << ",";
        statFile << stats.prevTotalTraversalCost << "," << stats.prevRelGap << ",";
        statFile << aonAssignment.stats.lastChecksum << std::endl;
      }

      if (verbose) {
        std::cout << "  Line search: " << stats.lastLineSearchTime << "ms";
        std::cout << "  Total: " << stats.lastRunningTime << "ms\n";
        std::cout << "  Prev total traversal cost: " << stats.prevTotalTraversalCost << "\n";
        std::cout << "  Prev relative gap: " << stats.prevRelGap << "\n";
        std::cout << std::flush;
      }
    }

    if (flowFile.is_open() && !outputIntermediates)
      FORALL_EDGES(graph, e) {
        const auto vol = trafficFlows[e];
        const auto sat = vol / graph.capacity(e);
        flowFile << aonAssignment.stats.numIterations << ',' << vol << ',' << sat << '\n';
      }

    if (distFile.is_open() && !outputIntermediates)
      for (const auto dist : aonAssignment.stats.lastDistances)
        distFile << aonAssignment.stats.numIterations << ',' << dist << '\n';

    if (verbose) {
      std::cout << "Total:\n";
      std::cout << "  Checksum: " << aonAssignment.stats.totalChecksum;
      std::cout << "  Prepro: " << aonAssignment.stats.totalPreprocessingTime << "ms";
      std::cout << "  Custom: " << aonAssignment.stats.totalCustomizationTime << "ms";
      std::cout << "  Queries: " << aonAssignment.stats.totalQueryTime << "ms";
      std::cout << "  Routing: " << aonAssignment.stats.totalRoutingTime << "ms\n";
      std::cout << "  Line search: " << stats.totalLineSearchTime << "ms";
      std::cout << "  Total: " << stats.totalRunningTime << "ms\n";
      std::cout << std::flush;
    }
  }

  // Returns the traffic flow on edge e.
  double trafficFlowOn(const int e) const {
    assert(e >= 0); assert(e < graph.numEdges());
    return trafficFlows[e];
  }

  FrankWolfeAssignmentStats stats; // Statistics about the execution.

 private:
  // Determines the initial solution.
  void determineInitialSolution(const int skipInterval) {
    #pragma omp parallel for schedule(static)
    FORALL_EDGES(graph, e)
      graph.traversalCost(e) = objFunction.derivative(e, 0);
    aonAssignment.run(skipInterval);
    #pragma omp parallel for schedule(static)
    FORALL_EDGES(graph, e)
      trafficFlows[e] = aonAssignment.trafficFlowOn(e);
  }

  // Updates traversal costs.
  void updateTraversalCosts() {
    auto totalTraversalCost = 0.0, totalPathCost = 0.0;
    #pragma omp parallel for reduction(+: totalTraversalCost, totalPathCost) schedule(static)
    FORALL_EDGES(graph, e) {
      graph.traversalCost(e) = objFunction.derivative(e, trafficFlows[e]);
      totalTraversalCost += trafficFlows[e] * traversalCostFunction(e, trafficFlows[e]);
      totalPathCost += trafficFlows[e] * graph.traversalCost(e);
    }
    stats.prevTotalTraversalCost = totalTraversalCost;
    stats.prevTotalPathCost = totalPathCost;
  }

  // Finds the descent direction.
  void findDescentDirection(const int skipInterval) {
    aonAssignment.run(skipInterval);
#ifndef TA_NO_CFW
    if (aonAssignment.stats.numIterations == 2) {
      FORALL_EDGES(graph, e)
        pointOfSight[e] = aonAssignment.trafficFlowOn(e);
      return;
    }

    auto num = 0.0, den = 0.0;
    #pragma omp parallel for reduction(+: num, den) schedule(static)
    FORALL_EDGES(graph, e) {
      const auto residualDirection = pointOfSight[e] - trafficFlows[e];
      const auto secondDerivative = objFunction.secondDerivative(e, trafficFlows[e]);
      const auto fwDirection = aonAssignment.trafficFlowOn(e) - trafficFlows[e];
      num += residualDirection * secondDerivative * fwDirection;
      den += residualDirection * secondDerivative * (fwDirection - residualDirection);
    }

    const auto alpha = std::min(std::max(0.0, num / den), 1 - 1e-15);
    #pragma omp parallel for schedule(static)
    FORALL_EDGES(graph, e)
      pointOfSight[e] = alpha * pointOfSight[e] + (1 - alpha) * aonAssignment.trafficFlowOn(e);
#endif
  }

  // Find the optimal move size.
  double findMoveSize() const {
    return bisectionMethod([this](const double tau) {
      auto sum = 0.0;
      #pragma omp parallel for reduction(+: sum) schedule(static)
      FORALL_EDGES(graph, e) {
#ifndef TA_NO_CFW
        const auto direction = pointOfSight[e] - trafficFlows[e];
#else
        const auto direction = aonAssignment.trafficFlowOn(e) - trafficFlows[e];
#endif
        sum += direction * objFunction.derivative(e, trafficFlows[e] + tau * direction);
      }
      return sum;
    }, 0, 1);
  }

  // Moves along the descent direction.
  void moveAlongDescentDirection(const double tau) {
    #pragma omp parallel for schedule(static)
    FORALL_EDGES(graph, e)
#ifndef TA_NO_CFW
      trafficFlows[e] += tau * (pointOfSight[e] - trafficFlows[e]);
#else
      trafficFlows[e] += tau * (aonAssignment.trafficFlowOn(e) - trafficFlows[e]);
#endif
  }

  using AonAssignment = AllOrNothingAssignment<ShortestPathAlgoT<Graph, TraversalCostAttribute>>;
  using TraversalCostFunction = TraversalCostFunctionT<Graph>;
  using ObjFunction = ObjFunctionT<TraversalCostFunction>;

  AonAssignment aonAssignment;                 // The all-or-nothing assignment subroutine.
  Graph& graph;                                // The network of interest.
  std::vector<double> trafficFlows;            // The traffic flows on the edges.
  std::vector<double> pointOfSight;            // The point defining the descent direction d = s - x
  TraversalCostFunction traversalCostFunction; // A functor returning the traversal cost of an edge.
  ObjFunction objFunction;                     // The objective function to be minimized (UE or SO).
  const bool verbose;                          // Should informative messages be displayed?
};

// An alias template for a user-equilibrium (UE) traffic assignment.
template <
    template <typename> class TraversalCostFunctionT,
    template <typename, typename> class ShortestPathAlgoT, typename GraphT>
using UEAssignment =
    FrankWolfeAssignment<UserEquilibrium, TraversalCostFunctionT, ShortestPathAlgoT, GraphT>;

// An alias template for a system-optimum (SO) traffic assignment.
template <
    template <typename> class TraversalCostFunctionT,
    template <typename, typename> class ShortestPathAlgoT, typename GraphT>
using SOAssignment =
    FrankWolfeAssignment<SystemOptimum, TraversalCostFunctionT, ShortestPathAlgoT, GraphT>;
