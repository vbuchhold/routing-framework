#pragma once

#include <algorithm>
#include <cassert>
#include <vector>

#include "Tools/Constants.h"

// Implementation of an addressable k-heap. It maintains a set of elements, each with an associated
// ID and key, under the standard priority queue operations. The elements are addressed by the IDs.
// This class is implemented as a min-heap, but can be easily turned into a max-heap by multiplying
// the keys by -1.
template <int K>
class AddressableKHeap {
  static_assert(K > 0, "parameter k must be strictly positive");

 public:
  // Constructs an empty addressable k-heap that can maintain elements with IDs from 0 to n - 1.
  explicit AddressableKHeap(const int n) {
    resize(n);
  }

  // Returns true if this heap contains no elements.
  bool empty() const {
    return heap.empty();
  }

  // Returns true if this heap contains an element with the specified ID.
  bool contains(const int id) const {
    assert(0 <= id); assert(id < elementIdToHeapIndex.size());
    return elementIdToHeapIndex[id] != INVALID_INDEX;
  }

  // Returns the ID of an element with minimum key.
  int minId() const {
    assert(!empty());
    return heap[0].id;
  }

  // Returns the minimum key.
  int minKey() const {
    assert(!empty());
    return heap[0].key;
  }

  // Removes all of the elements from this heap.
  void clear() {
    for (const auto& heapElement : heap)
      elementIdToHeapIndex[heapElement.id] = INVALID_INDEX;
    heap.clear();
  }

  // Ensures that this heap can maintain elements with IDs from 0 to n - 1.
  void resize(const int n) {
    clear();
    heap.reserve(n);
    elementIdToHeapIndex.assign(n, INVALID_INDEX);
  }

  // Inserts an element with the specified ID and key into this heap.
  void insert(const int id, const int key) {
    assert(!contains(id));
    heap.emplace_back(id, key);
    siftUp(heap.size() - 1);
  }

  // Returns the ID and key of an element with minimum key.
  void min(int& id, int& key) const {
    id = minId();
    key = minKey();
  }

  // Extracts an element with minimum key from this heap.
  void deleteMin(int& id, int& key) {
    assert(!empty());
    min(id, key);
    elementIdToHeapIndex[id] = INVALID_INDEX;
    heap.front() = heap.back();
    heap.pop_back();
    if (!empty())
      siftDown(0);
  }

  // Decreases the key of the element with the specified ID to newKey.
  void decreaseKey(const int id, const int newKey) {
    assert(contains(id));
    const int idx = elementIdToHeapIndex[id];
    assert(newKey <= heap[idx].key);
    heap[idx].key = newKey;
    siftUp(idx);
  }

 private:
  // An element in this heap, with an associated ID and key.
  struct HeapElement {
    // Constructs a heap element with the specified ID and key.
    HeapElement(const int id, const int key) : id(id), key(key) {}

    int id;
    int key;
  };

  // Moves the heap element stored in index idx toward the root until the heap property holds.
  void siftUp(int idx) {
    assert(0 <= idx); assert(idx < heap.size());
    const HeapElement elementToBeSiftedUp = heap[idx];
    while (idx > 0 && heap[getParent(idx)].key > elementToBeSiftedUp.key) {
      move(getParent(idx), idx);
      idx = getParent(idx);
    }
    heap[idx] = elementToBeSiftedUp;
    elementIdToHeapIndex[elementToBeSiftedUp.id] = idx;
  }

  // Moves the heap element stored in index idx down the tree until the heap property holds.
  void siftDown(int idx) {
    assert(0 <= idx); assert(idx < heap.size());
    const HeapElement elementToBeSiftedDown = heap[idx];
    while (getFirstChild(idx) < heap.size()) {
      const int firstChild = getFirstChild(idx);
      const int lastChild = std::min(getFirstChild(idx + 1), static_cast<int>(heap.size()));
      int minChild = firstChild;
      for (int i = firstChild + 1; i < lastChild; ++i)
        if (heap[i].key < heap[minChild].key)
          minChild = i;
      if (elementToBeSiftedDown.key > heap[minChild].key) {
        move(minChild, idx);
        idx = minChild;
      } else {
        break;
      }
    }
    heap[idx] = elementToBeSiftedDown;
    elementIdToHeapIndex[elementToBeSiftedDown.id] = idx;
  }

  // Returns the index of the parent of the specified child.
  static int getParent(const int idx) {
    return (idx - 1) / K;
  }

  // Returns the index of the first child of the specified parent.
  static int getFirstChild(const int idx) {
    return (idx * K) + 1;
  }

  // Moves the heap element stored in index idx1 to index idx2.
  void move(const int idx1, const int idx2) {
    assert(0 <= idx1); assert(idx1 < heap.size());
    assert(0 <= idx2); assert(idx2 < heap.size());
    heap[idx2] = heap[idx1];
    elementIdToHeapIndex[heap[idx2].id] = idx2;
  }

  std::vector<HeapElement> heap;         // A vector of all heap elements, being heap-ordered.
  std::vector<int> elementIdToHeapIndex; // A map from element IDs to heap indices.
};

// Aliases for several standard heaps.
using AddressableBinaryHeap = AddressableKHeap<2>;
using AddressableQuadheap = AddressableKHeap<4>;
