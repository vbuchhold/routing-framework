#pragma once

#include <array>
#include <cassert>
#include <vector>

#include "Algorithms/Dijkstra/BiDijkstra.h"
#include "Algorithms/Dijkstra/Dijkstra.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "DataStructures/Labels/SimdLabelSet.h"
#include "Tools/Constants.h"

namespace trafficassignment {

// An adapter that makes bidirectional search usable in the all-or-nothing assignment procedure.
template <typename InputGraphT, typename WeightT>
class BiDijkstraAdapter {
 public:
  // The label sets used by the standard and centralized bidirectional search.
  using LabelSet = BasicLabelSet<0, ParentInfo::FULL_PARENT_INFO>;
#if TA_LOG_K < 2 || defined(TA_NO_SIMD_SEARCH)
  using CentralizedLabelSet = BasicLabelSet<TA_LOG_K, ParentInfo::FULL_PARENT_INFO>;
#else
  using CentralizedLabelSet = SimdLabelSet<TA_LOG_K, ParentInfo::FULL_PARENT_INFO>;
#endif
  using InputGraph = InputGraphT;

  // The number of simultaneous shortest-path computations.
  static constexpr int K = CentralizedLabelSet::K;

  // Constructs an adapter for bidirectional search.
  BiDijkstraAdapter(const InputGraph& graph)
      : inputGraph(graph),
        reverseGraph(graph.getReverseGraph()),
        search(graph, reverseGraph),
        centralizedSearch(graph, reverseGraph) {
      assert(graph.isDefrag());
  }

  // Invoked before the first iteration.
  void preprocess() { /* do nothing */ }

  // Invoked before each iteration.
  void customize() {
    // Copy the current edge weights from the input graph to the reverse graph.
    FORALL_EDGES(reverseGraph, e) {
      const int weight = inputGraph.template get<WeightT>(reverseGraph.edgeId(e));
      reverseGraph.template get<WeightT>(e) = weight;
    }
  }

  // Computes the shortest path from s to t.
  void query(const int s, const int t) {
    search.run(s, t);
  }

  // Computes shortest paths from each source to its target simultaneously.
  void query(const std::array<int, K>& sources, const std::array<int, K>& targets) {
    centralizedSearch.run(sources, targets);
  }

  // Returns the length of the shortest path.
  int getDistance(const int /*dst*/) {
    return search.getDistance();
  }

  // Returns the length of the i-th centralized shortest path.
  int getDistance(const int /*dst*/, const int i) {
    return centralizedSearch.getDistance(i);
  }

  // Returns the edges on the (packed) shortest path.
  std::vector<int> getPackedEdgePath(const int /*dst*/) {
    return search.getEdgePath();
  }

  // Return the edges on the i-th (packed) centralized shortest path.
  std::vector<int> getPackedEdgePath(const int /*dst*/, const int i) {
    return centralizedSearch.getEdgePath(i);
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

  const InputGraph& inputGraph;                    // The input graph.
  InputGraph reverseGraph;                         // The reverse graph.
  BiDijkstra<Search> search;                       // Bidirectional search computing a single path.
  BiDijkstra<CentralizedSearch> centralizedSearch; // Bidirectional search computing multiple paths.
};

}
