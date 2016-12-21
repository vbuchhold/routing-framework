#pragma once

#include <algorithm>
#include <cassert>
#include <vector>

#include "Tools/CompilerSpecific.h"
#include "Tools/Constants.h"

// A container maintaining distance labels. It stores a global clock and a timestamp for each
// distance label. The timestamp indicates whether a distance label has a valid value or not.
template <typename DistanceLabelT>
class StampedDistanceLabelContainer {
 public:
  // Constructs a distance label container using timestamps.
  explicit StampedDistanceLabelContainer(const int numVertices)
      : distanceLabels(numVertices, INFTY), timestamps(numVertices, 0), clock(0) {}

  // Initializes all distance labels to infinity.
  void init() {
    ++clock;
    if (UNLIKELY(clock < 0)) {
      // Clock overflow occurred. Extremely unlikely.
      std::fill(distanceLabels.begin(), distanceLabels.end(), INFTY);
      std::fill(timestamps.begin(), timestamps.end(), 0);
      clock = 0;
    }
  }

  // Returns the distance label of v.
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
  std::vector<DistanceLabelT> distanceLabels; // The distance labels of the vertices.
  std::vector<int> timestamps;                // The timestamps indicating whether a label is valid.
  int clock;                                  // The global clock.
};
