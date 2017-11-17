#pragma once

#include <cassert>
#include <fstream>
#include <utility>
#include <vector>

#include <routingkit/contraction_hierarchy.h>

#include "DataStructures/Graph/Attributes/EdgeIdAttribute.h"
#include "DataStructures/Utilities/Permutation.h"
#include "Tools/BinaryIO.h"

// An implementation of a contraction hierarchy (CH).
template <typename SearchGraphT, typename WeightT>
class ContractionHierarchy {
  // The search graphs are required to have IDs associated with the edges.
  static_assert(SearchGraphT::template has<EdgeIdAttribute>(), "Search graph is missing edge IDs.");

  // Converts the specified CH from RoutingKit's representation into our representation.
  template <typename ContractionHierarchyT>
  friend ContractionHierarchyT convert(const RoutingKit::ContractionHierarchy&, const int);

 public:
  using SearchGraph = SearchGraphT; // The type of the upward and downward graph.
  using Weight = WeightT;           // The attribute used as edge weight.

  // Constructs an empty CH.
  ContractionHierarchy() : numOrigEdges(0) {}

  // Constructs a CH from a binary file.
  explicit ContractionHierarchy(std::ifstream& in) {
    readFrom(in);
  }

  // Returns the graph containing all upward edges.
  const SearchGraph& getUpwardGraph() const {
    return upwardGraph;
  }

  // Returns the graph containing all downward edges.
  const SearchGraph& getDownwardGraph() const {
    return downwardGraph;
  }

  // Returns the i-th vertex in the contraction order.
  int vertexInContractionOrder(const int i) const {
    assert(i >= 0); assert(i < order.size());
    return order[i];
  }

  // Returns the rank of vertex v.
  int rank(const int v) const {
    assert(v >= 0); assert(v < ranks.size());
    return ranks[v];
  }

  // Returns the number of shortcut edges.
  int numShortcuts() const {
    return constituentEdges.size();
  }

  // Returns true if edge e is a shortcut edge.
  bool isEdgeShortcut(const int e) const {
    return e >= numOrigEdges;
  }

  // Returns the first constituent edge of shortcut s.
  int shortcutsFirstEdge(const int s) const {
    assert(s >= numOrigEdges); assert(s < numOrigEdges + constituentEdges.size());
    return constituentEdges[s - numOrigEdges].first;
  }

  // Returns the second constituent edge of shortcut s.
  int shortcutsSecondEdge(const int s) const {
    assert(s >= numOrigEdges); assert(s < numOrigEdges + constituentEdges.size());
    return constituentEdges[s - numOrigEdges].second;
  }

  // Reads a CH from a binary file.
  void readFrom(std::ifstream& in) {
    upwardGraph.readFrom(in);
    downwardGraph.readFrom(in);
    order.readFrom(in);
    bio::read(in, numOrigEdges);
    bio::read(in, constituentEdges);
    ranks = order.getInversePermutation();
  }

  // Writes a CH to a binary file.
  void writeTo(std::ofstream& out) const {
    upwardGraph.writeTo(out);
    downwardGraph.writeTo(out);
    order.writeTo(out);
    bio::write(out, numOrigEdges);
    bio::write(out, constituentEdges);
  }

 private:
  SearchGraph upwardGraph;   // The graph containing all upward edges.
  SearchGraph downwardGraph; // The graph containing all downward edges.
  Permutation order;         // order[i] = v indicates that v was the i-th vertex contracted.
  Permutation ranks;         // ranks[v] = i indicates that v was the i-th vertex contracted.
  int numOrigEdges;          // The number of original edges in the input graph.

  // For each shortcut s, the two (shortcut) edges constituting s.
  std::vector<std::pair<int, int>> constituentEdges;
};
