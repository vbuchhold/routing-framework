#pragma once

#include <cassert>

#include "DataStructures/Graph/Attributes/AbstractAttribute.h"

// An attribute associating a travel time with each edge of a graph.
class TravelTimeAttribute : public AbstractAttribute<int> {
 public:
  static constexpr Type DEFAULT_VALUE = INFTY;         // The attribute's default value.
  static constexpr const char* NAME   = "travel_time"; // The attribute's unique name.

  // The travel times are stored internally as integers in units of 1/PRECISION seconds.
  static constexpr int PRECISION = 10;

  // Returns the travel time in 0.1 seconds of edge e.
  Type travelTime(const int e) const {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

  // Sets the travel time of edge e to the specified value in 0.1 seconds.
  void setTravelTime(const int e, const Type val) {
    assert(e >= 0); assert(e < values.size());
    assert(val >= 0);
    values[e] = val;
  }
};
