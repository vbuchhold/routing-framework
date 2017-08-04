#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <vector>

#include "DataStructures/Utilities/OriginDestination.h"
#include "Stats/TrafficAssignment/AllOrNothingAssignmentStats.h"
#include "Tools/CommandLine/ProgressBar.h"
#include "Tools/Simd/AlignVector.h"
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
  AllOrNothingAssignment(const InputGraph& graph,
                         const std::vector<ClusteredOriginDestination>& odPairs,
                         const bool verbose = true)
      : stats(odPairs.size()),
        shortestPathAlgo(graph),
        inputGraph(graph),
        odPairs(odPairs),
        verbose(verbose) {
    Timer timer;
    shortestPathAlgo.preprocess();
    stats.totalPreprocessingTime = timer.elapsed();
    stats.lastRoutingTime = stats.totalPreprocessingTime;
    stats.totalRoutingTime = stats.totalPreprocessingTime;
    if (verbose) std::cout << "  Prepro: " << stats.totalPreprocessingTime << "ms" << std::endl;
  }

  // Assigns all OD-flows to their currently shortest paths.
  void run() {
    Timer timer;
    ++stats.numIterations;
    if (verbose) std::cout << "Iteration " << stats.numIterations << ": " << std::flush;
    shortestPathAlgo.customize();
    stats.lastCustomizationTime = timer.elapsed();

    timer.restart();
    ProgressBar bar(odPairs.size(), verbose);
    trafficFlows.clear();
    trafficFlows.resize(inputGraph.numEdges() + shortestPathAlgo.getNumShortcuts());
    stats.startIteration();
    for (int i = 0, k = 1; i < odPairs.size(); i += k, k = 1) {
      if (i + 1 >= odPairs.size() || !odPairs[i + 1].hasSameZones(odPairs[i])) {
        // Run a single shortest-path computation.
        shortestPathAlgo.query(odPairs[i].origin, odPairs[i].destination);
        const int dist = shortestPathAlgo.getDistance(odPairs[i].destination);
        const std::vector<int>& path = shortestPathAlgo.getPackedEdgePath(odPairs[i].destination);
        processResult(i, dist, path);
      } else {
        // Run multiple shortest-path computations simultaneously.
        std::array<int, K> sources;
        std::array<int, K> targets;
        sources.fill(odPairs[i].origin);
        targets.fill(odPairs[i].destination);
        for (; k < K && i + k < odPairs.size() && odPairs[i + k].hasSameZones(odPairs[i]); ++k) {
          sources[k] = odPairs[i + k].origin;
          targets[k] = odPairs[i + k].destination;
        }
        shortestPathAlgo.query(sources, targets);
        for (int j = 0; j < k; ++j) {
          const int dst = odPairs[i].destination;
          const int dist = shortestPathAlgo.getDistance(dst, j);
          const std::vector<int>& path = shortestPathAlgo.getPackedEdgePath(dst, j);
          processResult(i + j, dist, path);
        }
      }
      bar += k;
    }

    // Propagate traffic flows from shortcut to original edges.
    const int maxShortcutId = inputGraph.numEdges() + shortestPathAlgo.getNumShortcuts() - 1;
    for (int s = maxShortcutId; s >= inputGraph.numEdges(); --s) {
      trafficFlows[shortestPathAlgo.getShortcutsFirstEdge(s)] += trafficFlows[s];
      trafficFlows[shortestPathAlgo.getShortcutsSecondEdge(s)] += trafficFlows[s];
    }
    stats.lastQueryTime = timer.elapsed();
    stats.avgChangeInDistances /= odPairs.size();
    stats.finishIteration();

    if (verbose) {
      std::cout << " done.\n";
      std::cout << "  Checksum: " << stats.lastChecksum;
      std::cout << "  Custom: " << stats.lastCustomizationTime << "ms";
      std::cout << "  Queries: " << stats.lastQueryTime << "ms";
      std::cout << "  Routing: " << stats.lastRoutingTime << "ms\n";
      std::cout << std::flush;
    }
  }

  // Returns the traffic flow on edge e.
  const int& trafficFlowOn(const int e) const {
    assert(e >= 0); assert(e < inputGraph.numEdges());
    return trafficFlows[e];
  }

  AllOrNothingAssignmentStats stats; // Statistics about the execution.

 private:
  // The maximum number of simultaneous shortest-path computations.
  static constexpr int K = ShortestPathAlgoT::K;

  // Processes the result of the shortest-path computation for the i-th OD-pair.
  void processResult(const int i, const int dist, const std::vector<int>& path) {
    // Maintain the avg and max change in the OD-distances between the last two iterations.
    const double change = 1.0 * std::abs(dist - stats.lastDistances[i]) / stats.lastDistances[i];
    stats.lastChecksum += dist;
    stats.lastDistances[i] = dist;
    stats.avgChangeInDistances += std::max(0.0, change);
    stats.maxChangeInDistances = std::max(stats.maxChangeInDistances, change);

    // Assign the OD-flow to each edge on the computed path.
    for (const auto e : path) {
      assert(e >= 0); assert(e < trafficFlows.size());
      ++trafficFlows[e];
    }
  }

  using ODPairs = std::vector<ClusteredOriginDestination>;

  ShortestPathAlgoT shortestPathAlgo; // Algo computing shortest paths between OD-pairs.
  const InputGraph& inputGraph;       // The input graph.
  const ODPairs& odPairs;             // The OD-pairs to be assigned onto the graph.
  AlignVector<int> trafficFlows;      // The traffic flows on the edges.
  const bool verbose;                 // Should informative messages be displayed?
};
