#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <utility>
#include <vector>

#include <routingkit/filter.h>
#include <routingkit/id_mapper.h>
#include <routingkit/inverse_vector.h>
#include <routingkit/nested_dissection.h>

#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Partitioning/SeparatorTree.h"
#include "Tools/Bitwise.h"

// Implementation of recursive bisection using the partitioning algorithm Inertial Flow. The input
// graph is recursively cut into two pieces, such that the number of cut edges is minimized. The
// balance parameter determines the minimum size of each of the two pieces.
class RecursiveBisection {
 public:
  // Returns a separator tree for the specified graph, computed by recursive bisection.
  template <typename GraphT>
  SeparatorTree run(const GraphT& graph, const double balance = 0.3) const {
    // Convert the graph to RoutingKit's graph representation.
    std::vector<float> lats(graph.numVertices());
    std::vector<float> lngs(graph.numVertices());
    std::vector<unsigned int> tails(graph.numEdges());
    std::vector<unsigned int> heads(graph.numEdges());
    FORALL_VERTICES(graph, u) {
      lats[u] = graph.latLng(u).latInDeg();
      lngs[u] = graph.latLng(u).lngInDeg();
      FORALL_INCIDENT_EDGES(graph, u, e) {
        tails[e] = u;
        heads[e] = graph.edgeHead(e);
      }
    }

    // Compute a separator tree by recursively cutting the graph into two pieces.
    SeparatorTree tree(graph.numVertices());
    auto subgraph = RoutingKit::make_graph_fragment(graph.numVertices(), tails, heads);
    recurseOnSubgraphs(subgraph, -1, lats, lngs, std::round(balance * 100), tree);
    return tree;
  }

 private:
  // Recursively bisects the connected components of the specified graph.
  void recurseOnSubgraphs(
      const RoutingKit::GraphFragment& graph, const int level,
      const std::vector<float>& lats, const std::vector<float>& lngs,
      const int balance, SeparatorTree& tree) const {
    // Decompose the graph into connected components and assign them to one of the two subgraphs.
    std::array<std::vector<RoutingKit::GraphFragment>, 2> subgraphs;
    for (auto& comp : decompose_graph_fragment_into_connected_components(graph)) {
      assert(comp.node_count() > 0);
      const bool belongsTo2ndSubgraph = tree.getSide(comp.global_node_id[0], std::max(level, 0));
      subgraphs[belongsTo2ndSubgraph].push_back(std::move(comp));
    }

    for (auto& subgraph : subgraphs) {
      const int msb = highestOneBit(subgraph.size() - 1);
      for (int i = 0; i < subgraph.size(); ++i) {
        auto& comp = subgraph[i];
        if (comp.node_count() == 1) {
          for (int j = 0; j < msb + 1; ++j)
            tree.setSide(comp.global_node_id[0], level + 1 + j, getBit(i, j));
        } else if (comp.node_count() == 2) {
          for (int v = 0; v < 2; ++v) {
            for (int j = 0; j < msb + 1; ++j)
              tree.setSide(comp.global_node_id[v], level + 1 + j, getBit(i, j));
            tree.setSide(comp.global_node_id[v], level + msb + 2, v);
          }
        } else {
          // Bisect the connected component.
          auto cut = inertial_flow(comp, balance, lats, lngs);
          RoutingKit::pick_smaller_side(cut);
          for (int v = 0; v < comp.node_count(); ++v) {
            for (int j = 0; j < msb + 1; ++j)
              tree.setSide(comp.global_node_id[v], level + 1 + j, getBit(i, j));
            tree.setSide(comp.global_node_id[v], level + msb + 2, cut.is_node_on_side.is_set(v));
          }

          // Remove all cut edges from the component, decomposing it into two (or more) subgraphs.
          const auto isNoCutEdge = RoutingKit::make_bit_vector(
              comp.arc_count(), [&comp, &cut](const int e) {
                assert(e >= 0); assert(e < comp.arc_count());
                const int tail = comp.tail[e];
                const int head = comp.head[e];
                return cut.is_node_on_side.is_set(tail) == cut.is_node_on_side.is_set(head);
              });
          inplace_keep_element_of_vector_if(isNoCutEdge, comp.tail);
          inplace_keep_element_of_vector_if(isNoCutEdge, comp.head);
          inplace_keep_element_of_vector_if(isNoCutEdge, comp.back_arc);
          comp.first_out = RoutingKit::invert_vector(comp.tail, comp.node_count());
          RoutingKit::LocalIDMapper map(isNoCutEdge);
          for (auto& e : comp.back_arc)
            e = map.to_local(e);

          recurseOnSubgraphs(comp, level + msb + 2, lats, lngs, balance, tree);
        }
      }
    }
  }
};
