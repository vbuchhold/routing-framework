#pragma once

#include <fstream>

#include "DataStructures/Utilities/Permutation.h"

// Implementation of a contraction hierarchy (CH).
template <typename GraphT, template <typename> class GetWeightT>
class ContractionHierarchy {
  // Only the CH preprocessing is allowed to construct new CHs.
  template <typename, template <typename> class>
  friend class CHPreprocessing;

 public:
  using Graph = GraphT;            // The type of the upward and downward graphs.
  template <typename G>
  using GetWeight = GetWeightT<G>; // A functor returning the edge weight used for routing.

  // Constructs a CH from a binary file.
  explicit ContractionHierarchy(std::ifstream& in) {
    readFrom(in);
  };

  // Returns the graph containing all upward edges.
  const Graph& upwardGraph() const {
    return upwardSearchGraph;
  }

  // Returns the graph containing all downward edges.
  const Graph& downwardGraph() const {
    return downwardSearchGraph;
  }

  // Returns the order in which vertices were contracted.
  const std::vector<int>& contractionOrder() const {
    return order;
  }

  // Reorders the vertices according to the specified permutation.
  void permuteVertices(const Permutation& perm) {
    upwardSearchGraph.permuteVertices(perm);
    downwardSearchGraph.permuteVertices(perm);
  }

  // Reads a CH from a binary file.
  void readFrom(std::ifstream& in) {
    upwardSearchGraph.readFrom(in);
    downwardSearchGraph.readFrom(in);
  }

  // Writes a CH to a binary file.
  void writeTo(std::ofstream& out) const {
    upwardSearchGraph.writeTo(out);
    downwardSearchGraph.writeTo(out);
  }

 private:
  // Constructs an empty CH.
  ContractionHierarchy() = default;

  Graph upwardSearchGraph;   // The graph containing all upward edges.
  Graph downwardSearchGraph; // The graph containing all downward edges.
  std::vector<int> order;    // order[i] = v indicates that v was the i-th vertex contracted.
};
