#pragma once

#include <cassert>
#include <cstdint>
#include <vector>

#include "Tools/Constants.h"

// This class represents a subset of a finite set of size n. Inserting elements, removing elements,
// and testing elements for membership take constant time. Iterating through and clearing a subset
// of size k both take time O(k).
class Subset {
 public:
  // Constructs an empty subset of a finite set of the specified size.
  explicit Subset(const int size) : elementsToIndices(size, INVALID_INDEX) {
    elements.reserve(size);
  }

  // Returns an iterator referring to the first element in the subset.
  std::vector<int32_t>::const_iterator begin() const noexcept {
    return elements.begin();
  }

  // Returns the past-the-end iterator for the subset.
  std::vector<int32_t>::const_iterator end() const noexcept {
    return elements.end();
  }

  // Returns the number of elements in the subset.
  int size() const noexcept {
	  return elements.size();
  }

  // Inserts the specified element into the subset. Invalidates only the past-the-end iterator.
  bool insert(const int element) {
    if (contains(element))
      return false;
    elementsToIndices[element] = elements.size();
    elements.push_back(element);
    return true;
  }

  // Removes the specified element from the subset. May invalidate all iterators.
  bool remove(const int element) {
    if (!contains(element))
      return false;
    elements[elementsToIndices[element]] = elements.back();
    elementsToIndices[elements.back()] = elementsToIndices[element];
    elements.pop_back();
    elementsToIndices[element] = INVALID_INDEX;
    return true;
  }

  // Removes all elements in the subset. May invalidate all iterators.
  void clear() {
    for (const auto element : elements)
      elementsToIndices[element] = INVALID_INDEX;
    elements.clear();
  }

  // Returns true if the subset contains the specified element.
  bool contains(const int element) const {
    assert(element >= 0); assert(element < elementsToIndices.size());
    return elementsToIndices[element] != INVALID_INDEX;
  }

 private:
  std::vector<int32_t> elements;          // The elements contained in the subset.
  std::vector<int32_t> elementsToIndices; // The index in the element array of each element.
};
