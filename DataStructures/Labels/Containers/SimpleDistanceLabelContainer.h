#pragma once

#include <algorithm>
#include <cassert>
#include <vector>

#include "Tools/Simd/AlignedVector.h"
#include "Tools/Constants.h"

// The simplest implementation of a container maintaining distance labels during a shortest-path
// search. The distance labels are explicitly initialized with a linear sweep over all labels.
template <typename DistanceLabelT>
class SimpleDistanceLabelContainer {
 public:
  // Constructs a distance label container with explicit initialization.
  explicit SimpleDistanceLabelContainer(const int numLabels) {
    resize(numLabels);
  }

  // Ensures that this container can hold the specified number of distance labels.
  void resize(const int numLabels) {
    const auto currentSize = distanceLabels.size();
    if (numLabels < currentSize)
      distanceLabels.erase(distanceLabels.begin() + numLabels, distanceLabels.end());
    else
      distanceLabels.insert(distanceLabels.end(), numLabels - currentSize, INFTY);
  }

  // Initializes all distance labels to infinity.
  void init() {
    std::fill(distanceLabels.begin(), distanceLabels.end(), INFTY);
  }

  // Returns a reference to the distance label of vertex v.
  DistanceLabelT& operator[](const int v) {
    assert(v >= 0); assert(v < distanceLabels.size());
    return distanceLabels[v];
  }

 private:
  AlignedVector<DistanceLabelT> distanceLabels; // The distance labels of the vertices.
};
