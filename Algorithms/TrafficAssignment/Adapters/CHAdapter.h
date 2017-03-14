#pragma once

#include <vector>

#include "Algorithms/CH/CHPreprocessing.h"
#include "Algorithms/CH/CHQuery.h"
#include "DataStructures/Graph/Attributes/EdgeIdAttribute.h"
#include "DataStructures/Graph/Attributes/RoutingCostAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "DataStructures/Utilities/OriginDestination.h"

// An adapter that makes CHs usable in the all-or-nothing assignment procedure.
template <typename InputGraphT, typename WeightT>
class CHAdapter {
 public:
  using InputGraph = InputGraphT;

  // Constructs an adapter for CHs.
  CHAdapter(const InputGraph& graph) : inputGraph(graph), chPrepro(graph), chSearch(ch) {}

  // Invoked before the first iteration.
  void preprocess() { /* do nothing */ }

  // Invoked before each iteration.
  void customize() {
    ch = chPrepro.run();
    chSearch.resize();
  }

  // Computes the shortest path between the specified OD-pair.
  void query(const OriginDestination& od) {
    chSearch.run(ch.rank(od.origin), ch.rank(od.destination));
  }

  // Returns the length of the shortest path computed last.
  int getDistance(const int /*dst*/) {
    return chSearch.getDistance();
  }

  // Returns the edges on the (packed) shortest path computed last.
  std::vector<int> getPackedEdgePath(const int /*dst*/) {
    return chSearch.getPackedEdgePath();
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
  using CHQuery = StandardCHQuery<CH, BasicLabelSet<1, ParentInfo::FULL_PARENT_INFO>>;

  const InputGraph& inputGraph; // The input graph.
  CH ch;                        // The CH rebuilt in each iteration.
  CHPrepro chPrepro;            // A CH preprocessing instance.
  CHQuery chSearch;             // A CH search instance.
};
