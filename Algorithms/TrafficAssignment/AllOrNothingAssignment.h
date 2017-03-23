#pragma once

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <vector>

#include "DataStructures/Utilities/OriginDestination.h"
#include "Stats/TrafficAssignment/AllOrNothingAssignmentStats.h"
#include "Tools/CommandLine/ProgressBar.h"
#include "Tools/Timer.h"

// Implementation of an iterative all-or-nothing traffic assignment. Each OD-pair is processed in
// turn and the corresponding OD-flow (in our case always a single flow unit) is assigned to each
// edge on the shortest path between O and D. Other O-D paths are not assigned any flow. The
// procedure can be used with different shortest-path algorithms.
template <typename ShortestPathAlgoT>
class AllOrNothingAssignment {
 private:
  using InputGraph = typename ShortestPathAlgoT::InputGraph;

 public:
  // Constructs an all-or-nothing assignment instance.
  AllOrNothingAssignment(const InputGraph& graph, const std::vector<OriginDestination>& odPairs,
                         const bool verbose = true)
      : stats(odPairs.size()),
        shortestPathAlgo(graph),
        inputGraph(graph),
        odPairs(odPairs),
        verbose(verbose) {
    Timer timer;
    shortestPathAlgo.preprocess();
    stats.totalPreprocessingTime = timer.elapsed();
    if (verbose) std::cout << "  Prepro: " << stats.totalPreprocessingTime << "ms" << std::endl;
  }

  // Assigns all OD-flows to their currently shortest paths.
  void run() {
    Timer timer;
    if (verbose) std::cout << "Iteration " << ++stats.numIterations << ": " << std::flush;
    shortestPathAlgo.customize();
    stats.lastCustomizationTime = timer.elapsed();

    timer.restart();
    ProgressBar bar(odPairs.size(), verbose);
    trafficFlows.clear();
    trafficFlows.resize(inputGraph.numEdges() + shortestPathAlgo.getNumShortcuts());
    stats.startIteration();
    for (int i = 0; i != odPairs.size(); ++i) {
      shortestPathAlgo.query(odPairs[i]);
      const int dist = shortestPathAlgo.getDistance(odPairs[i].destination);
      stats.lastChecksum += dist;

      // Keep track of the max and avg change in the OD-distances between the last two iterations.
      const double change = 1.0 * std::abs(dist - stats.lastDistances[i]) / stats.lastDistances[i];
      stats.lastDistances[i] = dist;
      stats.maxChangeInDistances = std::max(stats.maxChangeInDistances, change);
      stats.avgChangeInDistances += std::max(0.0, change);

      // Assign the OD-flow to each edge on the computed path.
      for (const auto e : shortestPathAlgo.getPackedEdgePath(odPairs[i].destination)) {
        assert(e >= 0); assert(e < trafficFlows.size());
        ++trafficFlows[e];
      }
      ++bar;
    }

    // Propagate traffic flows from shortcut to original edges.
    const int maxShortcutId = inputGraph.numEdges() + shortestPathAlgo.getNumShortcuts() - 1;
    for (int s = maxShortcutId; s >= inputGraph.numEdges(); --s) {
      trafficFlows[shortestPathAlgo.getShortcutsFirstEdge(s)] += trafficFlows[s];
      trafficFlows[shortestPathAlgo.getShortcutsSecondEdge(s)] += trafficFlows[s];
    }
    stats.lastQueryTime = timer.elapsed();
    stats.finishIteration();

    if (verbose) {
      std::cout << " done.\n";
      std::cout << "  Checksum: " << stats.lastChecksum;
      std::cout << "  Custom: " << stats.lastCustomizationTime << "ms";
      std::cout << "  Queries: " << stats.lastQueryTime << "ms";
      std::cout << "  Total: " << stats.lastCustomizationTime + stats.lastQueryTime << "ms\n";
      std::cout << std::flush;
    }
  }

  // Returns the traffic flow on edge e.
  int trafficFlowOn(const int e) const {
    assert(e >= 0); assert(e < inputGraph.numEdges());
    return trafficFlows[e];
  }

  AllOrNothingAssignmentStats stats; // Statistics about the execution.

 private:
  ShortestPathAlgoT shortestPathAlgo;            // Algo computing shortest paths between OD-pairs.
  const InputGraph& inputGraph;                  // The input graph.
  const std::vector<OriginDestination>& odPairs; // The OD-pairs to be assigned onto the graph.
  std::vector<int> trafficFlows;                 // The traffic flows on the edges.
  const bool verbose;                            // Should informative messages be displayed?
};
