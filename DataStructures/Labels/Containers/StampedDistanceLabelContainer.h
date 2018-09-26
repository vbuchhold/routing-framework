#pragma once

#include <algorithm>
#include <cassert>
#include <vector>

#include "Tools/Simd/AlignedVector.h"
#include "Tools/CompilerSpecific.h"
#include "Tools/Constants.h"

// A container maintaining distance labels. It stores a global clock and a timestamp for each
// distance label. The timestamp indicates whether a distance label has a valid value or not.
template <typename DistanceLabelT>
class StampedDistanceLabelContainer {
 public:
  // Constructs a distance label container using timestamps.
  explicit StampedDistanceLabelContainer(const int numVertices) : clock(0) {
    resize(numVertices);
  }

  // Ensures that this container can hold the specified number of distance labels.
  void resize(const int numVertices) {
    distanceLabels.resize(numVertices);
    timestamps.resize(numVertices);
  }

  // Initializes all distance labels to infinity.
  void init() {
    ++clock;
    if (UNLIKELY(clock < 0)) {
      // Clock overflow occurred. Extremely unlikely.
      std::fill(timestamps.begin(), timestamps.end(), 0);
      clock = 1;
    }
  }

  // Returns a reference to the distance label of v.
  DistanceLabelT& operator[](const int v) {
    assert(v >= 0); assert(v < distanceLabels.size());
    if (timestamps[v] != clock) {
      assert(timestamps[v] < clock);
      distanceLabels[v] = INFTY;
      timestamps[v] = clock;
    }
    return distanceLabels[v];
  }

 private:
  AlignedVector<DistanceLabelT> distanceLabels; // The distance labels of the vertices.
  std::vector<int> timestamps;                  // The timestamps indicating if a label is valid.
  int clock;                                    // The global clock.
};
