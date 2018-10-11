#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <vector>

#include <routingkit/nested_dissection.h>

#include "Algorithms/CCH/CCH.h"
#include "Algorithms/CCH/CCHMetric.h"
#include "Algorithms/CCH/EliminationTreeQuery.h"
#include "Algorithms/CH/CH.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "DataStructures/Labels/SimdLabelSet.h"
#include "DataStructures/Partitioning/SeparatorDecomposition.h"
#include "Tools/Simd/AlignedVector.h"

namespace trafficassignment {

// An adapter that makes CCHs usable in the all-or-nothing assignment procedure.
template <typename InputGraphT, typename WeightT>
class CCHAdapter {
 public:
#if TA_LOG_K < 2 || defined(TA_NO_SIMD_SEARCH)
  using LabelSet = BasicLabelSet<TA_LOG_K, ParentInfo::FULL_PARENT_INFO>;
#else
  using LabelSet = SimdLabelSet<TA_LOG_K, ParentInfo::FULL_PARENT_INFO>;
#endif
  using InputGraph = InputGraphT;

  static constexpr int K = LabelSet::K; // The number of simultaneous shortest-path computations.

  // The search algorithm using the graph and possibly auxiliary data to compute shortest paths.
  // Multiple instances can work on the same data concurrently.
  class QueryAlgo {
   public:
    // Constructs a query algorithm instance working on the specified data.
    QueryAlgo(
        const CH& minimumWeightedCH, const std::vector<int32_t>& eliminationTree,
        AlignedVector<int>& flowsOnUpEdges, AlignedVector<int>& flowsOnDownEdges)
        : minimumWeightedCH(minimumWeightedCH),
          search(minimumWeightedCH, eliminationTree),
          flowsOnUpEdges(flowsOnUpEdges),
          flowsOnDownEdges(flowsOnDownEdges),
          localFlowsOnUpEdges(flowsOnUpEdges.size()),
          localFlowsOnDownEdges(flowsOnDownEdges.size()) {
      assert(minimumWeightedCH.upwardGraph().numEdges() == flowsOnUpEdges.size());
      assert(minimumWeightedCH.downwardGraph().numEdges() == flowsOnDownEdges.size());
    }

    // Computes shortest paths from each source to its target simultaneously.
    void run(std::array<int, K>& sources, std::array<int, K>& targets, const int k) {
      // Run a centralized CH search.
      for (auto i = 0; i < K; ++i) {
        sources[i] = minimumWeightedCH.rank(sources[i]);
        targets[i] = minimumWeightedCH.rank(targets[i]);
      }
      search.run(sources, targets);

      // Assign flow to the edges on the computed paths.
      for (auto i = 0; i < k; ++i) {
        for (const auto e : search.getUpEdgePath(i)) {
          assert(e >= 0); assert(e < localFlowsOnUpEdges.size());
          ++localFlowsOnUpEdges[e];
        }
        for (const auto e : search.getDownEdgePath(i)) {
          assert(e >= 0); assert(e < localFlowsOnDownEdges.size());
          ++localFlowsOnDownEdges[e];
        }
      }
    }

    // Returns the length of the i-th shortest path.
    int getDistance(const int /*dst*/, const int i) {
      return search.getDistance(i);
    }

    // Adds the local flow counters to the global ones. Must be synchronized externally.
    void addLocalToGlobalFlows() {
      FORALL_EDGES(minimumWeightedCH.upwardGraph(), e)
        flowsOnUpEdges[e] += localFlowsOnUpEdges[e];
      FORALL_EDGES(minimumWeightedCH.downwardGraph(), e)
        flowsOnDownEdges[e] += localFlowsOnDownEdges[e];
    }

   private:
    const CH& minimumWeightedCH;            // The CH resulting from perfect customization.
    EliminationTreeQuery<LabelSet> search;  // The CH search on the minimum weighted CH.
    AlignedVector<int>& flowsOnUpEdges;     // The flows in the upward graph.
    AlignedVector<int>& flowsOnDownEdges;   // The flows in the downward graph.
    std::vector<int> localFlowsOnUpEdges;   // The local flows in the upward graph.
    std::vector<int> localFlowsOnDownEdges; // The local flows in the downward graph.
  };

  // Constructs an adapter for CCHs.
  explicit CCHAdapter(const InputGraph& inputGraph)
      : inputGraph(inputGraph), currentMetric(cch, &inputGraph.template get<WeightT>(0)) {
    assert(inputGraph.numEdges() > 0); assert(inputGraph.isDefrag());
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

    // Compute a separator decomposition for the input graph.
    const auto graph = RoutingKit::make_graph_fragment(inputGraph.numVertices(), tails, heads);
    auto computeSep = [&](const RoutingKit::GraphFragment& fragment) {
      const auto cut = inertial_flow(fragment, 30, lats, lngs);
      return derive_separator_from_cut(fragment, cut.is_node_on_side);
    };
    const auto decomp = compute_separator_decomposition(graph, computeSep);

    // Convert the separator decomposition to our representation.
    SeparatorDecomposition sepDecomp;
    for (const auto& n : decomp.tree) {
      SeparatorDecomposition::Node node;
      node.leftChild = n.left_child;
      node.rightSibling = n.right_sibling;
      node.firstSeparatorVertex = n.first_separator_vertex;
      node.lastSeparatorVertex = n.last_separator_vertex;
      sepDecomp.tree.push_back(node);
    }
    sepDecomp.order.assign(decomp.order.begin(), decomp.order.end());

    // Build the CCH.
    cch.preprocess(inputGraph, sepDecomp);
  }

  // Invoked before each iteration.
  void customize() {
    minimumWeightedCH = currentMetric.buildMinimumWeightedCH();
    flowsOnUpEdges.assign(minimumWeightedCH.upwardGraph().numEdges(), 0);
    flowsOnDownEdges.assign(minimumWeightedCH.downwardGraph().numEdges(), 0);
  }

  // Returns an instance of the query algorithm.
  QueryAlgo getQueryAlgoInstance() {
    return {minimumWeightedCH, cch.getEliminationTree(), flowsOnUpEdges, flowsOnDownEdges};
  }

  // Propagates the flows on the edges in the search graphs to the edges in the input graph.
  void propagateFlowsToInputEdges(AlignedVector<int>& flowsOnInputEdges) {
    const auto& upGraph = minimumWeightedCH.upwardGraph();
    const auto& downGraph = minimumWeightedCH.downwardGraph();
    for (auto u = inputGraph.numVertices() - 1; u >= 0; --u) {
      FORALL_INCIDENT_EDGES(upGraph, u, e)
        if (upGraph.unpackingInfo(e).second == INVALID_EDGE) {
          flowsOnInputEdges[upGraph.unpackingInfo(e).first] = flowsOnUpEdges[e];
        } else {
          flowsOnDownEdges[upGraph.unpackingInfo(e).first] += flowsOnUpEdges[e];
          flowsOnUpEdges[upGraph.unpackingInfo(e).second] += flowsOnUpEdges[e];
        }
      FORALL_INCIDENT_EDGES(downGraph, u, e)
        if (downGraph.unpackingInfo(e).second == INVALID_EDGE) {
          flowsOnInputEdges[downGraph.unpackingInfo(e).first] = flowsOnDownEdges[e];
        } else {
          flowsOnDownEdges[downGraph.unpackingInfo(e).first] += flowsOnDownEdges[e];
          flowsOnUpEdges[downGraph.unpackingInfo(e).second] += flowsOnDownEdges[e];
        }
    }
  }

 private:
  const InputGraph& inputGraph; // The input graph.
  CCH cch;                      // The metric-independent CCH.
  CCHMetric currentMetric;      // The current metric for the CCH.
  CH minimumWeightedCH;         // The minimum weighted CH resulting from perfect customization.

  AlignedVector<int> flowsOnUpEdges;   // The flows on the edges in the upward graph.
  AlignedVector<int> flowsOnDownEdges; // The flows on the edges in the downward graph.
};

}
