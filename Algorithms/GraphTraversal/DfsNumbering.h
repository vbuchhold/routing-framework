#pragma once

#include <vector>
#include <utility>

#include "Algorithms/GraphTraversal/DepthFirstSearch.h"

// An algorithm that numbers the vertices of a graph in the order in which they are reached during
// a DFS. The vertices may be reordered according to a DFS numbering in order to improve spatial
// locality during shortest-path computations. This class is an instantiation of the general DFS
// template, implementing the needed hook functions.
class DfsNumbering {
 public:
  // Constructs a DFS instance computing DFS numbers.
  DfsNumbering() : dfs(*this), verticesReached(0) {}

  // Returns a vector indexed by vertices with zero-based DFS numbers.
  template <typename GraphT>
  std::vector<int> run(const GraphT& graph) {
    dfs.run(graph);
    return std::vector<int>(std::move(dfsNumbers));
  }

  // Returns a vector indexed by vertices with zero-based DFS numbers.
  template <typename GraphT>
  std::vector<int> run(const GraphT& graph, const int s) {
    dfs.run(graph, s);
    return std::vector<int>(std::move(dfsNumbers));
  }

 private:
  // The DFS implementation is allowed to call the hook functions.
  friend class DepthFirstSearch<DfsNumbering>;

  // Functions for marking and unmarking vertices as reached.
  bool hasBeenReached(const int v) const { return dfsNumbers[v] >= 0; }
  void markAsReached(const int /*v*/)    { /* implicitly marked as reached by hook functions */ }
  void unmarkVertices(const int numVertices) { dfsNumbers.assign(numVertices, -1); }

  // Hook functions called during the execution of the DFS.
  void init()                                         { verticesReached = 0; }
  void root(const int s)                              { dfsNumbers[s] = verticesReached++; }
  void traverseTreeEdge(const int /*v*/, const int w) { dfsNumbers[w] = verticesReached++; }
  void traverseNonTreeEdge(const int /*v*/, const int /*w*/) { /* not needed */ }
  void backtrack(const int /*u*/, const int /*v*/)           { /* not needed */ }

  DepthFirstSearch<DfsNumbering> dfs; // The DFS instance executing the actual search.
  std::vector<int> dfsNumbers;        // The DFS numbers of the vertices.
  int verticesReached;                // The number of vertices reached so far.
};
