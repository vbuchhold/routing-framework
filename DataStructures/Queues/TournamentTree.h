#pragma once

#include <array>
#include <cassert>
#include <utility>

#include "Tools/Math.h"

// A tournament tree, also known as selection tree or loser tree. Suppose we have k sorted sequences
// that are to be merged into a single output sequence. We repeatedly have to find the smallest from
// the leading elements in the k sequences. This can be done efficiently by a tournament tree.
template <int logK>
class TournamentTree {
 public:
  static constexpr int K = 1 << logK; // The number of sequences that are to be merged.

  // Constructs a tournament tree given the leading elements in the k sequences.
  explicit TournamentTree(const std::array<int, K>& keys) {
    build(keys);
  }

  // Builds a tournament tree given the leading elements in the k sequences.
  void build(const std::array<int, K>& keys) {
    std::array<Vertex, K> winners;
    for (int i = K - 2; i >= 0; i -= 2)
      minmax({keys[i], i}, {keys[i + 1], i + 1}, winners[getParent(i + K)], tree[getParent(i + K)]);
    for (int i = K - 2; i > 0; i -= 2)
      minmax(winners[i], winners[i + 1], winners[getParent(i)], tree[getParent(i)]);
    tree[0] = K != 1 ? winners[1] : Vertex{keys[0], 0};
  }

  // Returns the key of the smallest element.
  int minKey() const {
    return tree[0].key;
  }

  // Returns the index of the sequence where the smallest element comes from.
  int minSeq() const {
    return tree[0].seq;
  }

  // Deletes the smallest element and inserts the next element from the corresponding sequence.
  void deleteMin(const int keyOfNextElement) {
    Vertex& winner = tree[0];
    winner.key = keyOfNextElement;
    int parent = winner.seq + K;
    for (int i = 0; i < logK; ++i) {
      parent = getParent(parent);
      if (tree[parent] < winner)
        std::swap(tree[parent], winner);
    }
  }

 private:
  // A vertex in the tournament tree. It represents the leading element from one of the sequences.
  struct Vertex {
    // Returns true if lhs wins against rhs.
    friend bool operator<(const Vertex& lhs, const Vertex& rhs) {
      return lhs.key < rhs.key;
    }

    int key; // The key of the element.
    int seq; // The index of the sequence where the element comes from.
  };

  // Returns the index of the parent of the specified child.
  static constexpr int getParent(const int child) {
    assert(child > 0); assert(child < 2 * K);
    return child / 2;
  }

  std::array<Vertex, K> tree; // The actual tournament tree.
};
