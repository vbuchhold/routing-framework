#pragma once

#include <cassert>
#include <cstdint>
#include <fstream>
#include <vector>

#include "DataStructures/Utilities/Permutation.h"
#include "Tools/BinaryIO.h"

// A separator decomposition of an undirected n-vertex graph G = (V, E) is a rooted tree T whose
// nodes are disjoint subsets of V and that is recursively defined as follows. If n = 1, then T
// consists of a single node equal to V. If n > 1, then T consists of a root whose removal separates
// G into several subgraphs. The children of the root are the roots of separator decompositions of
// these subgraphs.
struct SeparatorDecomposition {
  // A node in the separator decomposition.
  struct Node {
    int32_t leftChild;            // The index of the left child.
    int32_t rightSibling;         // The index of the right sibling.
    int32_t firstSeparatorVertex; // The index of the first separator vertex.
    int32_t lastSeparatorVertex;  // The index one past the last separator vertex.
  };

  // Returns the index in the tree of the left child of this node.
  int leftChild(const int node) const noexcept {
    assert(node >= 0); assert(node < tree.size());
    return tree[node].leftChild;
  }

  // Returns the index in the tree of the right sibling of this node.
  int rightSibling(const int node) const noexcept {
    assert(node >= 0); assert(node < tree.size());
    return tree[node].rightSibling;
  }

  // Returns the index in the order of the first separator vertex contained in this node.
  int firstSeparatorVertex(const int node) const noexcept {
    assert(node >= 0); assert(node < tree.size());
    return tree[node].firstSeparatorVertex;
  }

  // Returns the index in the order one past the last separator vertex contained in this node.
  int lastSeparatorVertex(const int node) const noexcept {
    assert(node >= 0); assert(node < tree.size());
    return tree[node].lastSeparatorVertex;
  }

  // Reads the separator decomposition from the specified binary file.
  void readFrom(std::ifstream& in) {
    bio::read(in, tree);
    order.readFrom(in);
  }

  // Writes the separator decomposition to the specified binary file.
  void writeTo(std::ofstream& out) const {
    bio::write(out, tree);
    order.writeTo(out);
  }

  std::vector<Node> tree; // The rooted tree representing this separator decomposition.
  Permutation order;      // The associated nested dissection order.
};
