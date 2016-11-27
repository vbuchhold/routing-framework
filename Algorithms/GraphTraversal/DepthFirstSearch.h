#pragma once

#include <cassert>
#include <stack>
#include <vector>

#include "DataStructures/Graph/Graph.h"

// Non-recursive implementation of depth-first search. This class is an algorithm template.
// Concrete DFS-based algorithms (like computing a DFS numbering or strongly connected components)
// are derived from it by implementing various hook functions called during execution of the DFS.
template <typename ConcreteAlgoT>
class DepthFirstSearch {
  // Only concrete DFS-based algorithms are allowed to access the members of this class.
  friend ConcreteAlgoT;

 private:
  // An active vertex, i.e., a vertex that has been reached but not finished.
  struct ActiveVertex {
    // Constructs an active vertex.
    ActiveVertex(const int id, const int nextUnexploredEdge)
        : id(id), nextUnexploredEdge(nextUnexploredEdge) {}

    int id;                 // The ID of this active vertex.
    int nextUnexploredEdge; // The next unexplored incident edge.
  };

  // Constructs a DFS instance.
  explicit DepthFirstSearch(ConcreteAlgoT& algo) : concreteAlgo(algo) {
    // Push a sentinel onto the stack.
    activeVertices.emplace(-1, -1);
  }

  // Runs a non-recursive DFS.
  template <typename GraphT>
  void run(const GraphT& graph) {
    concreteAlgo.unmarkVertices(graph.numVertices());
    concreteAlgo.init();
    FORALL_VERTICES(graph, s)
      if (!concreteAlgo.hasBeenReached(s))
        growDfsTree(graph, s);
  }

  // Runs a non-recursive DFS starting at s.
  template <typename GraphT>
  void run(const GraphT& graph, const int s) {
    concreteAlgo.unmarkVertices(graph.numVertices());
    concreteAlgo.init();
    growDfsTree(graph, s);
  }

  // Grows a DFS tree rooted at s, calling hook functions for DFS-based algorithms during execution.
  template <typename GraphT>
  void growDfsTree(const GraphT& graph, const int s) {
    concreteAlgo.markAsReached(s);
    concreteAlgo.root(s);
    assert(activeVertices.size() == 1);
    activeVertices.top().id = s; // Set the sentinel.
    activeVertices.emplace(s, graph.firstEdge(s));

    // Explore all active vertices in the stack, except for the sentinel.
    while (activeVertices.size() > 1) {
      ActiveVertex &v = activeVertices.top();
      if (v.nextUnexploredEdge != graph.lastEdge(v.id)) {
        // There are edges out of v that have not been explored yet. Explore one of them.
        const int head = graph.edgeHead(v.nextUnexploredEdge);
        ++v.nextUnexploredEdge;
        if (concreteAlgo.hasBeenReached(head)) {
          concreteAlgo.traverseNonTreeEdge(v.id, head);
        } else {
          concreteAlgo.traverseTreeEdge(v.id, head);
          concreteAlgo.markAsReached(head);
          activeVertices.emplace(head, graph.firstEdge(head));
        }
      } else {
        // All edges out of v have been explored. Return from v along the incoming tree edge.
        const int head = v.id;
        activeVertices.pop();
        assert(activeVertices.size() == 1 || graph.containsEdge(activeVertices.top().id, head));
        concreteAlgo.backtrack(activeVertices.top().id, head);
      }
    }
  }

  using Stack = std::stack<ActiveVertex, std::vector<ActiveVertex>>;

  ConcreteAlgoT& concreteAlgo; // The concrete DFS-based algorithm.
  Stack activeVertices;        // The stack of active vertices.
};
