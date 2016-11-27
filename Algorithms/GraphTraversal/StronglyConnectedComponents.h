#pragma once

#include <cassert>
#include <stack>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "Algorithms/GraphTraversal/DepthFirstSearch.h"
#include "Tools/Workarounds.h"

// An algorithm that computes the strongly connected components of a graph. It implements the
// Cheriyan-Mehlhorn algorithm, also known as Gabow's algorithm, using optimization techniques from
// Mehlhorn et al. (2007). This class is an instantiation of the general DFS template, implementing
// the needed hook functions.
class StronglyConnectedComponents {
 public:
  // Constructs a DFS instance computing the SCCs of a graph.
  StronglyConnectedComponents()
      : dfs(*this), currentNegatedDfsNumber(-2), largestScc(-1), largestSccSize(0) {}

  // Computes the SCCs of the specified graph.
  template <typename GraphT>
  void run(const GraphT graph) {
    dfs.run(graph);
    for (const auto representative : components) {
      unused(representative);
      assert(representative >= 0); assert(representative < graph.numVertices());
      assert(components[representative] == representative);
    }
    assert(openRepresentatives.empty());
    assert(openVertices.empty());
    assert(currentNegatedDfsNumber == -graph.numVertices() - 2);
    assert(largestSccSize >= 0); assert(largestSccSize <= graph.numVertices());
  }

  // Returns a vector indexed by vertices where vec[u] = vec[v] iff u and v belong to the same SCC.
  std::vector<int> getStronglyConnectedComponents() const {
    return components;
  }

  // Returns a bitmask indexed by vertices where bitmask[v] is set iff v lies in the largest SCC.
  boost::dynamic_bitset<> getLargestSccAsBitmask() const {
    boost::dynamic_bitset<> bitmask(components.size());
    for (int i = 0; i < components.size(); ++i)
      bitmask[i] = components[i] == largestScc;
    return bitmask;
  }

 private:
  // The DFS implementation is allowed to call the hook functions.
  friend class DepthFirstSearch<StronglyConnectedComponents>;

  using StackT = std::stack<int, std::vector<int>>; // A stack of integers.

  // Functions for marking and unmarking vertices as reached.
  bool hasBeenReached(const int v) const { return components[v] != -1; }
  void markAsReached(const int /*v*/)    { /* implicitly marked as reached by hook functions */ }
  void unmarkVertices(const int numVertices) {
    // During the execution of the algorithm, the data member components stores three kinds of
    // information. Initially, components[v] stores -1, indicating that v has not yet been reached
    // by the DFS. When v is reached, components[v] is set to the negated DFS number of v. Finally,
    // when the SCC of v is closed, components[v] is set to the representative vertex of the SCC.
    // To disambiguate the values 0 and -1, DFS numbers start at 2.
    components.assign(numVertices, -1);
    currentNegatedDfsNumber = -2;
  }

  // Hook functions called during the execution of the DFS.
  void init()            { largestScc = -1; largestSccSize = 0; }
  void root(const int s) { traverseTreeEdge(s, s); }
  void traverseTreeEdge(const int /*v*/, const int w) {
    components[w] = currentNegatedDfsNumber--;
    // w becomes an open SCC on its own.
    openRepresentatives.push(components[w]);
    openVertices.push(w);
  }
  void traverseNonTreeEdge(const int /*v*/, const int w) {
    if (components[w] < -1)
      // w belongs to an open SCC. Merge all SCCs lying on the newly formed cycle into a single SCC.
      while (components[w] > openRepresentatives.top())
        openRepresentatives.pop();
  }
  void backtrack(const int /*u*/, const int v) {
    if (components[v] == openRepresentatives.top()) {
      // v is a representative vertex. Close its SCC.
      openRepresentatives.pop();
      int w, size = 0;
      do {
        w = openVertices.top();
        openVertices.pop();
        components[w] = v;
        ++size;
      } while (w != v);

      // Keep track of the largest SCC seen so far.
      if (size > largestSccSize) {
        largestScc = v;
        largestSccSize = size;
      }
    }
  }

  DepthFirstSearch<StronglyConnectedComponents> dfs; // DFS instance executing the actual search.

  std::vector<int> components; // Stores different kinds of information; see member unmarkVertices.
  StackT openRepresentatives;  // Stores the (DFS numbers of the) representatives of the open SCCs.
  StackT openVertices;         // Stores all vertices that belong to an open SCC.

  int currentNegatedDfsNumber; // Maintains the current (negated) DFS number during the DFS.
  int largestScc;              // The representative vertex of the largest SCC seen so far.
  int largestSccSize;          // The size of the largest SCC seen so far.
};
