#pragma once

#include <cassert>
#include <vector>

// This class maintains a partition of the set 0..n under the following operations: testing whether
// two elements are in the same subset and joining two subsets into one. This implementation uses
// three optimizations: union by rank, path compression, and storing parent pointers and ranks in
// a single array.
class UnionFind {
 public:
  // Constructs a partition with each element being a block on its own.
  explicit UnionFind(const int numElements) : parent(numElements, numElements) {}

  // Returns the representative of the block containing i.
  int find(int i) {
    assert(i >= 0); assert(i < parent.size());
    int rep = i;
    // Find the representative of i.
    while (parent[rep] < parent.size())
      rep = parent[rep];
    // Redirect the parent pointers.
    while (parent[i] < parent.size()) {
      const int tmp = parent[i];
      parent[i] = rep;
      i = tmp;
    }
    return rep;
  }

  // Joins the blocks containing i and j (i and j must be representatives of different blocks).
  void link(const int i, const int j) {
    assert(find(i) == i);
    assert(find(j) == j);
    assert(i != j);
    if (parent[i] < parent[j]) {
      parent[i] = j;
    } else {
      if (parent[i] == parent[j])
        ++parent[i];
      parent[j] = i;
    }
  }

  // Joins the blocks containing i and j.
  void unite(const int i, const int j) {
    if (find(i) != find(j))
      link(find(i), find(j));
  }

 private:
  std::vector<int> parent; // parent[i] stores parent for each nonrep i, and rank for each rep i.
};
