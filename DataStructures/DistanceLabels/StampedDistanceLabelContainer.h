#pragma once

#include <algorithm>
#include <cassert>
#include <vector>

#include "Tools/CompilerSpecific.h"
#include "Tools/Constants.h"

// A container maintaining distance labels. It stores a global clock and a timestamp for each
// distance label. The timestamp indicates whether a distance label has a valid value or not.
class StampedDistanceLabelContainer {
 public:
  // Constructs a distance label container using timestamps.
  explicit StampedDistanceLabelContainer(const int numVertices)
      : distanceLabels(numVertices, {INFTY, 0}), clock(0) {}

  // Initializes all distance labels to infinity.
  void init() {
    ++clock;
    if (clock < 0) {
      // Clock overflow occurred. Extremely unlikely.
      std::fill(distanceLabels.begin(), distanceLabels.end(), DistanceLabel{INFTY, 0});
      clock = 0;
    }
  }

  // Returns the tentative distance of v.
  int& operator[](const int v) {
    assert(v >= 0); assert(v < distanceLabels.size());
    if (distanceLabels[v].timestamp != clock) {
      assert(distanceLabels[v].timestamp < clock);
      distanceLabels[v].distance = INFTY;
      distanceLabels[v].timestamp = clock;
    }
    return distanceLabels[v].distance;
  }

 private:
  // A distance label for a particular vertex, maintaining the tentative distance and timestamp.
  struct DistanceLabel {
    int distance;
    int timestamp;
  };

  std::vector<DistanceLabel> distanceLabels; // The distance labels of the vertices.
  int clock;                                 // The global clock.
};
