#pragma once

#include <cassert>
#include <type_traits>
#include <utility>
#include <vector>

#include <routingkit/contraction_hierarchy.h>

#include "Algorithms/CH/ContractionHierarchy.h"
#include "DataStructures/Graph/Attributes/EdgeIdAttribute.h"
#include "DataStructures/Graph/Graph.h"

// Implementation of a CH preprocessing. Vertices are reordered in increasing rank order.
template <typename GraphT, template <typename> class GetWeightT>
class CHPreprocessing {
 public:
  using CH = ContractionHierarchy<GraphT, GetWeightT>; // The CH type we build.

  // Constructs a CH preprocessing instance.
  CHPreprocessing(const GraphT& inputGraph) : inputGraph(inputGraph) {
    assert(inputGraph.isDefrag());
  }

  // Contracts the input graph.
  CH run() {
    // Convert the input graph into RoutingKit's graph representation.
    std::vector<unsigned int> tails(inputGraph.numEdges());
    std::vector<unsigned int> heads(inputGraph.numEdges());
    std::vector<unsigned int> weights(inputGraph.numEdges());
    FORALL_VALID_EDGES(inputGraph, u, e) {
      tails[e] = u;
      heads[e] = inputGraph.edgeHead(e);
      weights[e] = getWeight(inputGraph, e);
    }

    const RoutingKit::ContractionHierarchy ch =
       RoutingKit::ContractionHierarchy::build(inputGraph.numVertices(), tails, heads, weights);

    // Convert the search graphs into our graph representation.
    GraphT upwardGraph = getSearchGraph(ch.forward);
    GraphT downwardGraph = getSearchGraph(ch.backward);

    setEdgeIds(upwardGraph, downwardGraph, ch);
    CH ret;
    ret.upwardSearchGraph = std::move(upwardGraph);
    ret.downwardSearchGraph = std::move(downwardGraph);
    ret.order = {ch.order.begin(), ch.order.end()};
    return ret;
  }

  // Returns the search graph of the specified side.
  GraphT getSearchGraph(const RoutingKit::ContractionHierarchy::Side& side) const {
    GraphT graph(inputGraph.numVertices(), side.head.size());
    for (int v = 0; v != inputGraph.numVertices(); ++v) {
      graph.appendVertex();
      for (int e = side.first_out[v]; e != side.first_out[v + 1]; ++e) {
        graph.appendEdge(side.head[e]);
        getWeight(graph, e) = side.weight[e];
      }
    }
    return graph;
  }

  // Assigns IDs to original and shortcut edges. The IDs of the original edges are given by their
  // indices in the original graph. Shortcuts have sequential IDs starting from m - the number of
  // of original edges -, given by the order in which they were added during preprocessing.

  template <typename G>
  std::enable_if_t<G::template has<EdgeIdAttribute>()>
  setEdgeIds(G& upwardGraph, G& downwardGraph, const RoutingKit::ContractionHierarchy& ch) const {
    int nextShortcutId = inputGraph.numEdges();
    FORALL_VERTICES(inputGraph, u) {
      FORALL_INCIDENT_EDGES(upwardGraph, u, e) {
        const bool isOrig = ch.forward.is_shortcut_an_original_arc.is_set(e);
        upwardGraph.edgeId(e) = isOrig ? ch.forward.shortcut_first_arc[e] : nextShortcutId++;
      }
      FORALL_INCIDENT_EDGES(downwardGraph, u, e) {
        const bool isOrig = ch.backward.is_shortcut_an_original_arc.is_set(e);
        downwardGraph.edgeId(e) = isOrig ? ch.backward.shortcut_first_arc[e] : nextShortcutId++;
      }
    }
  }

  template <typename G>
  std::enable_if_t<!G::template has<EdgeIdAttribute>()>
  setEdgeIds(G&, G&, const RoutingKit::ContractionHierarchy&) const {}

 private:
  const GraphT& inputGraph;     // The input graph that should be contracted.
  GetWeightT<GraphT> getWeight; // A functor returning the edge weight used for routing.
};
