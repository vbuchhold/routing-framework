#pragma once

#include <cassert>
#include <fstream>
#include <utility>

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
