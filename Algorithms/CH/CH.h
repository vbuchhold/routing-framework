#pragma once

#include <cassert>
#include <fstream>
#include <utility>

#include <routingkit/contraction_hierarchy.h>

#include "DataStructures/Graph/Attributes/TraversalCostAttribute.h"
#include "DataStructures/Graph/Attributes/UnpackingInfoAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Utilities/Permutation.h"

// A weighted contraction hierarchy. The contraction order is determined online and bottom-up.
class CH {
 public:
  // The type of the upward and downward graph.
  using EdgeAttributes = EdgeAttrs<TraversalCostAttribute, UnpackingInfoAttribute>;
  using SearchGraph = StaticGraph<VertexAttrs<>, EdgeAttributes>;

  // Constructs an empty CH.
  CH() = default;

  // Constructs a CH from the specified upward and downward graph.
  template <typename UpGraphT, typename DownGraphT, typename OrderT, typename RanksT>
  CH(UpGraphT&& upGraph, DownGraphT&& downGraph, OrderT&& order, RanksT&& ranks) noexcept
      : upGraph(std::forward<UpGraphT>(upGraph)),
        downGraph(std::forward<DownGraphT>(downGraph)),
        order(std::forward<OrderT>(order)),
        ranks(std::forward<RanksT>(ranks)) {
    assert(this->upGraph.numVertices() == this->downGraph.numVertices());
    assert(this->upGraph.numVertices() == this->order.size());
    assert(this->ranks == this->order.getInversePermutation());
  }

  // Constructs a CH from the specified binary file.
  explicit CH(std::ifstream& in) {
    readFrom(in);
  }

  // Builds a weighted CH for the specified graph.
  template <typename WeightT, typename InputGraphT>
  void preprocess(const InputGraphT& inputGraph) {
    const auto numVertices = inputGraph.numVertices();
    const auto numEdges = inputGraph.numEdges();
    std::vector<unsigned int> tails(numEdges);
    std::vector<unsigned int> heads(numEdges);
    std::vector<unsigned int> weights(numEdges);
    FORALL_VALID_EDGES(inputGraph, u, e) {
      tails[e] = u;
      heads[e] = inputGraph.edgeHead(e);
      weights[e] = inputGraph.template get<WeightT>(e);
    }
    const auto ch = RoutingKit::ContractionHierarchy::build(numVertices, tails, heads, weights);

    upGraph.clear();
    downGraph.clear();
    upGraph.reserve(numVertices, ch.forward.head.size());
    downGraph.reserve(numVertices, ch.backward.head.size());
    for (int u = 0; u < numVertices; ++u) {
      upGraph.appendVertex();
      downGraph.appendVertex();
      for (int e = ch.forward.first_out[u]; e < ch.forward.first_out[u + 1]; ++e) {
        upGraph.appendEdge(ch.forward.head[e]);
        upGraph.traversalCost(e) = ch.forward.weight[e];
        upGraph.unpackingInfo(e).first = ch.forward.shortcut_first_arc[e];
        upGraph.unpackingInfo(e).second = ch.forward.shortcut_second_arc[e];
        if (ch.forward.is_shortcut_an_original_arc.is_set(e))
          upGraph.unpackingInfo(e).second = INVALID_EDGE;
      }
      for (int e = ch.backward.first_out[u]; e < ch.backward.first_out[u + 1]; ++e) {
        downGraph.appendEdge(ch.backward.head[e]);
        downGraph.traversalCost(e) = ch.backward.weight[e];
        downGraph.unpackingInfo(e).first = ch.backward.shortcut_first_arc[e];
        downGraph.unpackingInfo(e).second = ch.backward.shortcut_second_arc[e];
        if (ch.backward.is_shortcut_an_original_arc.is_set(e))
          downGraph.unpackingInfo(e).second = INVALID_EDGE;
      }
    }

    order.assign(ch.order.begin(), ch.order.end());
    ranks.assign(ch.rank.begin(), ch.rank.end());
  }

  // Returns the upward graph.
  const SearchGraph& upwardGraph() const noexcept {
    return upGraph;
  }

  // Returns the downward graph.
  const SearchGraph& downwardGraph() const noexcept {
    return downGraph;
  }

  // Returns the i-th vertex in the contraction order.
  int contractionOrder(const int i) const noexcept {
    assert(i >= 0); assert(i < order.size());
    return order[i];
  }

  // Returns the position of vertex v in the contraction order
  int rank(const int v) const noexcept {
    assert(v >= 0); assert(v < ranks.size());
    return ranks[v];
  }

  // Reads the CH from the specified binary file.
  void readFrom(std::ifstream& in) {
    upGraph.readFrom(in);
    downGraph.readFrom(in);
    order.readFrom(in);
    ranks.readFrom(in);
  }

  // Writes the CH to the specified binary file.
  void writeTo(std::ofstream& out) const {
    upGraph.writeTo(out);
    downGraph.writeTo(out);
    order.writeTo(out);
    ranks.writeTo(out);
  }

 private:
  SearchGraph upGraph;   // The upward graph.
  SearchGraph downGraph; // The downward graph.
  Permutation order;     // order[i] = v indicates that v was the i-th vertex contracted.
  Permutation ranks;     // ranks[v] = i indicates that v was the i-th vertex contracted.
};
