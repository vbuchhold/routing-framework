#pragma once

#include <cassert>
#include <vector>

#include "Algorithms/Dijkstra/Dijkstra.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "DataStructures/Utilities/OriginDestination.h"
#include "Tools/Constants.h"

namespace trafficassignment {

// An adapter that makes Dijkstra's algorithm usable in the all-or-nothing assignment procedure.
template <typename InputGraphT, typename WeightT>
class DijkstraAdapter {
 public:
  using InputGraph = InputGraphT;

  // Constructs an adapter for Dijkstra's algorithm.
  DijkstraAdapter(const InputGraph& graph) : dijkstra(graph) {}

  // Invoked before the first iteration.
  void preprocess() { /* do nothing */ }

  // Invoked before each iteration.
  void customize() { /* do nothing */ }

  // Computes the shortest path between the specified OD-pair.
  void query(const OriginDestination& od) {
    dijkstra.run(od.origin, od.destination);
  }

  // Returns the length of the shortest path computed last.
  int getDistance(const int dst) {
    return dijkstra.getDistance(dst);
  }

  // Returns the edges on the (packed) shortest path computed last.
  std::vector<int> getPackedEdgePath(const int dst) {
    return dijkstra.getReverseEdgePath(dst);
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
  using LabelSet = BasicLabelSet<0, ParentInfo::FULL_PARENT_INFO>;

  StandardDijkstra<InputGraph, WeightT, LabelSet> dijkstra; // The Dijkstra search.
};

}
