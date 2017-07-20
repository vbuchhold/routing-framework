#pragma once

#include <array>
#include <vector>

#include "Algorithms/CH/CHPreprocessing.h"
#include "Algorithms/CH/CHQuery.h"
#include "DataStructures/Graph/Attributes/EdgeIdAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "DataStructures/Labels/SimdLabelSet.h"
#include "Tools/Constants.h"

namespace trafficassignment {

// An adapter that makes CHs usable in the all-or-nothing assignment procedure.
template <typename InputGraphT, typename WeightT>
class CHAdapter {
 public:
  // The label sets used by the standard and centralized CH search.
  using LabelSet = BasicLabelSet<0, ParentInfo::FULL_PARENT_INFO>;
#if TA_LOG_K < 2 || defined(TA_NO_SIMD_SEARCH)
  using CentralizedLabelSet = BasicLabelSet<TA_LOG_K, ParentInfo::FULL_PARENT_INFO>;
#else
  using CentralizedLabelSet = SimdLabelSet<TA_LOG_K, ParentInfo::FULL_PARENT_INFO>;
#endif
  using InputGraph = InputGraphT;

  // The number of simultaneous shortest-path computations.
  static constexpr int K = CentralizedLabelSet::K;

  // Constructs an adapter for CHs.
  CHAdapter(const InputGraph& graph)
      : inputGraph(graph),
        chPrepro(graph),
        search(ch),
        centralizedSearch(ch) {}

  // Invoked before the first iteration.
  void preprocess() { /* do nothing */ }

  // Invoked before each iteration.
  void customize() {
    ch = chPrepro.run();
    search.resize();
    centralizedSearch.resize();
  }

  // Computes the shortest path from s to t.
  void query(const int s, const int t) {
    search.run(ch.rank(s), ch.rank(t));
  }

  // Computes shortest paths from each source to its target simultaneously.
  void query(std::array<int, K>& sources, std::array<int, K>& targets) {
    for (int i = 0; i < K; ++i) {
      sources[i] = ch.rank(sources[i]);
      targets[i] = ch.rank(targets[i]);
    }
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
    return search.getPackedEdgePath();
  }

  // Return the edges on the i-th (packed) centralized shortest path.
  std::vector<int> getPackedEdgePath(const int /*dst*/, const int i) {
    return centralizedSearch.getPackedEdgePath(i);
  }

  // Returns the first constituent edge of shortcut s.
  int getShortcutsFirstEdge(const int s) const {
    return ch.shortcutsFirstEdge(s);
  }

  // Returns the second constituent edge of shortcut s.
  int getShortcutsSecondEdge(const int s) const {
    return ch.shortcutsSecondEdge(s);
  }

  // Returns the number of shortcut edges.
  int getNumShortcuts() const {
    return ch.numShortcuts();
  }

 private:
  using CHGraph = StaticGraph<VertexAttrs<>, EdgeAttrs<EdgeIdAttribute, WeightT>>;
  using CHPrepro = CHPreprocessing<InputGraph, CHGraph, WeightT>;
  using CH = typename CHPrepro::CH;
  using Search = StandardCHQuery<CH, LabelSet>;
  using CentralizedSearch = StandardCHQuery<CH, CentralizedLabelSet>;

  const InputGraph& inputGraph;        // The input graph.
  CH ch;                               // The CH rebuilt in each iteration.
  CHPrepro chPrepro;                   // A CH preprocessing instance.
  Search search;                       // CH search computing a single path.
  CentralizedSearch centralizedSearch; // CH search computing multiple paths.
};

}
