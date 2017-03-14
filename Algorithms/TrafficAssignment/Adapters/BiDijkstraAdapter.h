#pragma once

#include <cassert>
#include <vector>

#include "Algorithms/Dijkstra/BiDijkstra.h"
#include "Algorithms/Dijkstra/Dijkstra.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "DataStructures/Utilities/OriginDestination.h"
#include "Tools/Constants.h"

// An adapter that makes bidirectional search usable in the all-or-nothing assignment procedure.
template <typename InputGraphT, typename WeightT>
class BiDijkstraAdapter {
 public:
  using InputGraph = InputGraphT;

  // Constructs an adapter for Dijkstra's algorithm.
  BiDijkstraAdapter(const InputGraph& graph)
      : inputGraph(graph),
        reverseGraph(graph.getReverseGraph()),
        biDijkstra(graph, reverseGraph) {
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

  // Computes the shortest path between the specified OD-pair.
  void query(const OriginDestination& od) {
    biDijkstra.run(od.origin, od.destination);
  }

  // Returns the length of the shortest path computed last.
  int getDistance(const int /*dst*/) {
    return biDijkstra.getDistance();
  }

  // Returns the edges on the (packed) shortest path computed last.
  std::vector<int> getPackedEdgePath(const int /*dst*/) {
    return biDijkstra.getEdgePath();
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
  using LabelSet = BasicLabelSet<1, ParentInfo::FULL_PARENT_INFO>;
  using Dijkstra = StandardDijkstra<InputGraph, WeightT, LabelSet>;

  const InputGraph& inputGraph;     // The input graph.
  InputGraph reverseGraph;          // The reverse graph.
  BiDijkstra<Dijkstra> biDijkstra;  // The bidirectional search.
};
