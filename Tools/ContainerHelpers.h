#pragma once

#include <algorithm>

// Returns true if the specified sequence container contains the specified element.
template <typename InputIteratorT, typename T>
inline bool contains(InputIteratorT first, InputIteratorT last, const T& val) {
  return std::find(first, last, val) != last;
}
