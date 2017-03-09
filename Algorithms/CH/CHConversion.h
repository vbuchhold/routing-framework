#pragma once

#include <routingkit/contraction_hierarchy.h>

#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Utilities/Permutation.h"

// Converts the specified CH from RoutingKit's representation to our representation.
template <typename ContractionHierarchyT>
ContractionHierarchyT convert(const RoutingKit::ContractionHierarchy& src, const int numOrigEdges) {
  using SearchGraph = typename ContractionHierarchyT::SearchGraph;
  typename ContractionHierarchyT::template GetWeight<SearchGraph> getWeight;

  ContractionHierarchyT dst;
  dst.order = Permutation({src.order.begin(), src.order.end()});
  dst.ranks = Permutation({src.rank.begin(), src.rank.end()});
  dst.numOrigEdges = numOrigEdges;

  // Convert the upward graph to our graph representation.
  dst.upwardGraph.reserve(src.node_count(), src.forward.head.size());
  for (int v = 0; v != src.node_count(); ++v) {
    dst.upwardGraph.appendVertex();
    for (int e = src.forward.first_out[v]; e != src.forward.first_out[v + 1]; ++e) {
      dst.upwardGraph.appendEdge(src.forward.head[e]);
      getWeight(dst.upwardGraph, e) = src.forward.weight[e];
    }
  }

  // Convert the downward graph to our graph representation.
  dst.downwardGraph.reserve(src.node_count(), src.backward.head.size());
  for (int v = 0; v != src.node_count(); ++v) {
    dst.downwardGraph.appendVertex();
    for (int e = src.backward.first_out[v]; e != src.backward.first_out[v + 1]; ++e) {
      dst.downwardGraph.appendEdge(src.backward.head[e]);
      getWeight(dst.downwardGraph, e) = src.backward.weight[e];
    }
  }

  // Assign IDs to original and shortcut edges. The IDs of the original edges are given by their
  // indices in the original graph. Shortcuts have sequential IDs starting from m (the number of
  // original edges), given by the order in which they were added during preprocessing.
  int nextShortcutId = numOrigEdges;
  FORALL_VERTICES(dst.upwardGraph, u) {
    FORALL_INCIDENT_EDGES(dst.upwardGraph, u, e) {
      if (src.forward.is_shortcut_an_original_arc.is_set(e)) {
        // e is an original edge.
        dst.upwardGraph.edgeId(e) = src.forward.shortcut_first_arc[e];
      } else {
        // e is a shortcut.
        dst.upwardGraph.edgeId(e) = nextShortcutId++;
        const int firstEdge = dst.downwardGraph.edgeId(src.forward.shortcut_first_arc[e]);
        const int secondEdge = dst.upwardGraph.edgeId(src.forward.shortcut_second_arc[e]);
        dst.constituentEdges.emplace_back(firstEdge, secondEdge);
      }
    }
    FORALL_INCIDENT_EDGES(dst.downwardGraph, u, e) {
      if (src.backward.is_shortcut_an_original_arc.is_set(e)) {
        // e is an original edge.
        dst.downwardGraph.edgeId(e) = src.backward.shortcut_first_arc[e];
      } else {
        // e is a shortcut.
        dst.downwardGraph.edgeId(e) = nextShortcutId++;
        const int firstEdge = dst.downwardGraph.edgeId(src.backward.shortcut_first_arc[e]);
        const int secondEdge = dst.upwardGraph.edgeId(src.backward.shortcut_second_arc[e]);
        dst.constituentEdges.emplace_back(firstEdge, secondEdge);
      }
    }
  }

  return dst;
}
