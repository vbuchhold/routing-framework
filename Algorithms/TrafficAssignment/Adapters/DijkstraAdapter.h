#pragma once

#include <array>
#include <cassert>
#include <vector>

#include "Algorithms/Dijkstra/Dijkstra.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "DataStructures/Labels/SimdLabelSet.h"
#include "Tools/Constants.h"

namespace trafficassignment {

// An adapter that makes Dijkstra's algorithm usable in the all-or-nothing assignment procedure.
template <typename InputGraphT, typename WeightT>
class DijkstraAdapter {
 public:
  // The label sets used by the standard and centralized Dijkstra search.
  using LabelSet = BasicLabelSet<0, ParentInfo::FULL_PARENT_INFO>;
#if TA_LOG_K < 2 || defined(TA_NO_SIMD_SEARCH)
  using CentralizedLabelSet = BasicLabelSet<TA_LOG_K, ParentInfo::FULL_PARENT_INFO>;
#else
  using CentralizedLabelSet = SimdLabelSet<TA_LOG_K, ParentInfo::FULL_PARENT_INFO>;
#endif
  using InputGraph = InputGraphT;

  // The number of simultaneous shortest-path computations.
  static constexpr int K = CentralizedLabelSet::K;

  // Constructs an adapter for Dijkstra's algorithm.
  DijkstraAdapter(const InputGraph& graph) : search(graph), centralizedSearch(graph) {}

  // Invoked before the first iteration.
  void preprocess() { /* do nothing */ }

  // Invoked before each iteration.
  void customize() { /* do nothing */ }

  // Computes the shortest path from s to t.
  void query(const int s, const int t) {
    search.run(s, t);
  }

  // Computes shortest paths from each source to its target simultaneously.
  void query(const std::array<int, K>& sources, const std::array<int, K>& targets) {
    centralizedSearch.run(sources, targets);
  }

  // Returns the length of the shortest path.
  int getDistance(const int dst) {
    return search.getDistance(dst);
  }

  // Returns the length of the i-th centralized shortest path.
  int getDistance(const int dst, const int i) {
    return centralizedSearch.getDistance(dst, i);
  }

  // Returns the edges on the (packed) shortest path.
  std::vector<int> getPackedEdgePath(const int dst) {
    return search.getReverseEdgePath(dst);
  }

  // Return the edges on the i-th (packed) centralized shortest path.
  std::vector<int> getPackedEdgePath(const int dst, const int i) {
    return centralizedSearch.getReverseEdgePath(dst, i);
  }

  // Returns the first constituent edge of shortcut s.
  int getShortcutsFirstEdge(const int /*s*/) const {
    assert(false);
    return INVALID_EDGE;
  }

  // Returns the second constituent edge of shortcut s.
  int getShortcutsSecondEdge(const int /*s*/) const {
    assert(false);
    return INVALID_EDGE;
  }

  // Returns the number of shortcut edges.
  int getNumShortcuts() const {
    return 0;
  }

 private:
  using Search = StandardDijkstra<InputGraph, WeightT, LabelSet>;
  using CentralizedSearch = StandardDijkstra<InputGraph, WeightT, CentralizedLabelSet>;

  Search search;                       // Dijkstra search computing a single path.
  CentralizedSearch centralizedSearch; // Dijkstra search computing multiple paths.
};

}
