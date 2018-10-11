#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <vector>

#include "Algorithms/Dijkstra/BiDijkstra.h"
#include "Algorithms/Dijkstra/Dijkstra.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "DataStructures/Labels/SimdLabelSet.h"
#include "Tools/Simd/AlignedVector.h"

namespace trafficassignment {

// An adapter that makes bidirectional search usable in the all-or-nothing assignment procedure.
template <typename InputGraphT, typename WeightT>
class BiDijkstraAdapter {
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
        const InputGraph& inputGraph, const InputGraph& reverseGraph,
        AlignedVector<int>& flowsOnForwardEdges, AlignedVector<int>& flowsOnReverseEdges)
        : search(inputGraph, reverseGraph),
          flowsOnForwardEdges(flowsOnForwardEdges),
          flowsOnReverseEdges(flowsOnReverseEdges),
          localFlowsOnForwardEdges(flowsOnForwardEdges.size()),
          localFlowsOnReverseEdges(flowsOnReverseEdges.size()) {
      assert(inputGraph.numEdges() == flowsOnForwardEdges.size());
      assert(reverseGraph.numEdges() == flowsOnReverseEdges.size());
    }

    // Computes shortest paths from each source to its target simultaneously.
    void run(std::array<int, K>& sources, std::array<int, K>& targets, const int k) {
      // Run a centralized bidirectional search.
      search.run(sources, targets);

      // Assign flow to the edges on the computed paths.
      for (auto i = 0; i < k; ++i) {
        for (const auto e : search.getEdgePathToMeetingVertex(i)) {
          assert(e >= 0); assert(e < localFlowsOnForwardEdges.size());
          ++localFlowsOnForwardEdges[e];
        }
        for (const auto e : search.getEdgePathFromMeetingVertex(i)) {
          assert(e >= 0); assert(e < localFlowsOnReverseEdges.size());
          ++localFlowsOnReverseEdges[e];
        }
      }
    }

    // Returns the length of the i-th shortest path.
    int getDistance(const int /*dst*/, const int i) {
      return search.getDistance(i);
    }

    // Adds the local flow counters to the global ones. Must be synchronized externally.
    void addLocalToGlobalFlows() {
      for (auto e = 0; e < flowsOnForwardEdges.size(); ++e) {
        flowsOnForwardEdges[e] += localFlowsOnForwardEdges[e];
        flowsOnReverseEdges[e] += localFlowsOnReverseEdges[e];
      }
    }

   private:
    using Dijkstra = StandardDijkstra<InputGraph, WeightT, LabelSet>;

    BiDijkstra<Dijkstra> search;               // The bidirectional search.
    AlignedVector<int>& flowsOnForwardEdges;   // The flows in the forward graph.
    AlignedVector<int>& flowsOnReverseEdges;   // The flows in the reverse graph.
    std::vector<int> localFlowsOnForwardEdges; // The local flows in the forward graph.
    std::vector<int> localFlowsOnReverseEdges; // The local flows in the reverse graph.
  };

  // Constructs an adapter for bidirectional search.
  explicit BiDijkstraAdapter(const InputGraph& inputGraph)
      : inputGraph(inputGraph),
        reverseGraph(inputGraph.getReverseGraph()),
        flowsOnForwardEdges(inputGraph.numEdges()),
        flowsOnReverseEdges(inputGraph.numEdges()) {
    assert(inputGraph.isDefrag());
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
    std::fill(flowsOnForwardEdges.begin(), flowsOnForwardEdges.end(), 0);
    std::fill(flowsOnReverseEdges.begin(), flowsOnReverseEdges.end(), 0);
  }

  // Returns an instance of the query algorithm.
  QueryAlgo getQueryAlgoInstance() {
    return {inputGraph, reverseGraph, flowsOnForwardEdges, flowsOnReverseEdges};
  }

  // Propagates the flows on the edges in the search graphs to the edges in the input graph.
  void propagateFlowsToInputEdges(AlignedVector<int>& flowsOnInputEdges) {
    assert(flowsOnInputEdges.size() == flowsOnInputEdges.size());
    flowsOnInputEdges.swap(flowsOnForwardEdges);
    FORALL_EDGES(inputGraph, e)
      flowsOnInputEdges[reverseGraph.edgeId(e)] += flowsOnReverseEdges[e];
  }

 private:
  const InputGraph& inputGraph; // The input graph.
  InputGraph reverseGraph;      // The reverse graph.

  AlignedVector<int> flowsOnForwardEdges; // The flows on the edges in the forward graph.
  AlignedVector<int> flowsOnReverseEdges; // The flows on the edges in the reverse graph.
};

}
