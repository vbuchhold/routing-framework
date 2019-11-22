#pragma once

#include "DataStructures/Containers/BitVector.h"

// The class represents a subset of a finite set of size n. It is implemented as a bit vector of
// size n. Each bit indicates whether an element of the finite set is in the subset or not.
class Subset {
 public:
  // Constructs an (empty or complete) subset of a finite set of the specified size.
  explicit Subset(const int size, const bool complete = false)
      : isElementContained(size, complete) {}

  // Returns true if the subset contains the element with the specified index.
  bool contains(const int idx) const {
    return isElementContained[idx];
  }

  // Inserts the element with the specified index into the subset.
  void insert(const int idx) {
    isElementContained[idx] = true;
  }

  // Removes the element with the specified index idx from the subset.
  void remove(const int idx) {
    isElementContained[idx] = false;
  }

 private:
  BitVector isElementContained; // Each bit indicates whether an element is in the subset or not.
};
