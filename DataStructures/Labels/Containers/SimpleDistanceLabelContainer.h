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
    distanceLabels.resize(numLabels, INFTY);
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
