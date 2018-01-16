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

  // Some aliases concerning CHs.
  using CHGraph = StaticGraph<VertexAttrs<>, EdgeAttrs<EdgeIdAttribute, WeightT>>;
  using CHPrepro = CHPreprocessing<InputGraph, CHGraph, WeightT>;
  using CH = typename CHPrepro::CH;

  // The number of simultaneous shortest-path computations.
  static constexpr int K = CentralizedLabelSet::K;

  // The search algorithm using the graph and possibly auxiliary data to compute shortest paths.
  // Multiple instances may work on the same data concurrently.
  class QueryAlgo {
   public:
    // Constructs a query algorithm instance working on the specified data.
    explicit QueryAlgo(const CH& ch) : search(ch), centralizedSearch(ch), ch(ch) {}

    // Computes the shortest path from s to t.
    void run(const int s, const int t) {
      search.run(ch.rank(s), ch.rank(t));
    }

    // Computes shortest paths from each source to its target simultaneously.
    void run(std::array<int, K>& sources, std::array<int, K>& targets) {
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

    // Returns the edges on the i-th (packed) centralized shortest path.
    std::vector<int> getPackedEdgePath(const int /*dst*/, const int i) {
      return centralizedSearch.getPackedEdgePath(i);
    }

   private:
    using Search = StandardCHQuery<CH, LabelSet>;
    using CentralizedSearch = StandardCHQuery<CH, CentralizedLabelSet>;

    Search search;                       // CH search for a single path.
    CentralizedSearch centralizedSearch; // CH search for multiple paths.
    const CH& ch;                        // The CH rebuilt in each iteration.
  };

  // Constructs an adapter for CHs.
  CHAdapter(const InputGraph& graph) : chPrepro(graph) {}

  // Invoked before the first iteration.
  void preprocess() { /* do nothing */ }

  // Invoked before each iteration.
  void customize() {
    ch = chPrepro.run();
  }

  // Returns an instance of the query algorithm.
  QueryAlgo getQueryAlgoInstance() const {
    return QueryAlgo(ch);
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
  CHPrepro chPrepro; // A CH preprocessing instance.
  CH ch;             // The CH rebuilt in each iteration.
};

}
