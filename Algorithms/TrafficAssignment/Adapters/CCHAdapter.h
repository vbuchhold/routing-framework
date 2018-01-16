#pragma once

#include <array>
#include <cassert>
#include <vector>

#include <routingkit/customizable_contraction_hierarchy.h>
#include <routingkit/nested_dissection.h>

#include "Algorithms/CCH/EliminationTreeQuery.h"
#include "Algorithms/CH/CHConversion.h"
#include "Algorithms/CH/ContractionHierarchy.h"
#include "DataStructures/Graph/Attributes/EdgeIdAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "DataStructures/Labels/SimdLabelSet.h"
#include "Tools/Constants.h"

namespace trafficassignment {

// An adapter that makes CCHs usable in the all-or-nothing assignment procedure.
template <typename InputGraphT, typename WeightT>
class CCHAdapter {
 private:
  // The type of the CH resulting from perfectly customizing the CCH.
  using CHGraph = StaticGraph<VertexAttrs<>, EdgeAttrs<EdgeIdAttribute, WeightT>>;
  using CH = ContractionHierarchy<CHGraph, WeightT>;

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

  // The search algorithm using the graph and possibly auxiliary data to compute shortest paths.
  // Multiple instances may work on the same data concurrently.
  class QueryAlgo {
   public:
    // Constructs a query algorithm instance working on the specified data.
    QueryAlgo(const CH& perfectCH, const std::vector<int>& eliminationTree)
        : search(perfectCH, eliminationTree),
          centralizedSearch(perfectCH, eliminationTree),
          perfectCH(perfectCH) {}

    // Computes the shortest path from s to t.
    void run(const int s, const int t) {
      search.run(perfectCH.rank(s), perfectCH.rank(t));
    }

    // Computes shortest paths from each source to its target simultaneously.
    void run(std::array<int, K>& sources, std::array<int, K>& targets) {
      for (int i = 0; i < K; ++i) {
        sources[i] = perfectCH.rank(sources[i]);
        targets[i] = perfectCH.rank(targets[i]);
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

   private:
    using Search = EliminationTreeQuery<CH, LabelSet>;
    using CentralizedSearch = EliminationTreeQuery<CH, CentralizedLabelSet>;

    Search search;                       // CH search on the perfect CH for a single path.
    CentralizedSearch centralizedSearch; // CH search on the perfect CH for multiple paths.
    const CH& perfectCH;                 // The CH resulting from perfectly customizing the CCH.
  };

  // Constructs an adapter for CCHs.
  explicit CCHAdapter(const InputGraph& graph) : inputGraph(graph) {
    assert(graph.numEdges() > 0); assert(graph.isDefrag());
  }

  // Invoked before the first iteration.
  void preprocess() {
    // Convert the input graph to RoutingKit's graph representation.
    std::vector<float> lats(inputGraph.numVertices());
    std::vector<float> lngs(inputGraph.numVertices());
    std::vector<unsigned int> tails(inputGraph.numEdges());
    std::vector<unsigned int> heads(inputGraph.numEdges());
    FORALL_VERTICES(inputGraph, u) {
      lats[u] = inputGraph.latLng(u).latInDeg();
      lngs[u] = inputGraph.latLng(u).lngInDeg();
      FORALL_INCIDENT_EDGES(inputGraph, u, e) {
        tails[e] = u;
        heads[e] = inputGraph.edgeHead(e);
      }
    }

    // Compute a nested dissection order for the input graph.
    const auto graph = RoutingKit::make_graph_fragment(inputGraph.numVertices(), tails, heads);
    auto computeSeparator = [&lats, &lngs](const RoutingKit::GraphFragment& fragment) {
      const auto cut = inertial_flow(fragment, 30, lats, lngs);
      return derive_separator_from_cut(fragment, cut.is_node_on_side);
    };
    const auto order = compute_nested_node_dissection_order(graph, computeSeparator);

    // Build the metric-independent CCH.
    cch = RoutingKit::CustomizableContractionHierarchy(order, tails, heads);
    eliminationTree.assign(cch.elimination_tree_parent.begin(), cch.elimination_tree_parent.end());
    eliminationTree.back() = INVALID_VERTEX;
    const int* const weights = &inputGraph.template get<WeightT>(0);
    currentMetric = {cch, reinterpret_cast<const unsigned int*>(weights)};
  }

  // Invoked before each iteration.
  void customize() {
    // Customize the CCH using perfect witness searches.
    const int numEdges = inputGraph.numEdges();
    perfectCH = convert<CH>(
        currentMetric.build_contraction_hierarchy_using_perfect_witness_search(), numEdges);
  }

  // Returns an instance of the query algorithm.
  QueryAlgo getQueryAlgoInstance() const {
    return QueryAlgo(perfectCH, eliminationTree);
  }

  // Returns the first constituent edge of shortcut s.
  int getShortcutsFirstEdge(const int s) const {
    return perfectCH.shortcutsFirstEdge(s);
  }

  // Returns the second constituent edge of shortcut s.
  int getShortcutsSecondEdge(const int s) const {
    return perfectCH.shortcutsSecondEdge(s);
  }

  // Returns the number of shortcut edges.
  int getNumShortcuts() const {
    return perfectCH.numShortcuts();
  }

 private:
  using CCH = RoutingKit::CustomizableContractionHierarchy;
  using CCHMetric = RoutingKit::CustomizableContractionHierarchyMetric;

  const InputGraph& inputGraph;     // The input graph.
  CCH cch;                          // The metric-independent CCH.
  std::vector<int> eliminationTree; // eliminationTree[v] is the parent of v in the tree.
  CCHMetric currentMetric;          // The current metric for the CCH.
  CH perfectCH;                     // The CH resulting from perfectly customizing the CCH.
};

}
