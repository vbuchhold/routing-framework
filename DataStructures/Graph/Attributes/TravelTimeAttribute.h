#pragma once

#include <cassert>

#include "DataStructures/Graph/Attributes/AbstractAttribute.h"
#include "Tools/Constants.h"

// An attribute associating a travel time with each edge of a graph.
class TravelTimeAttribute : public AbstractAttribute<int> {
 public:
  // Returns the attribute's default value.
  static Type defaultValue() {
    return INFTY;
  }

  // Returns the travel time in tenths of a second on edge e.
  const Type& travelTime(const int e) const {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

  // Returns a reference to the travel time in tenths of a second on edge e.
  Type& travelTime(const int e) {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

 protected:
  static constexpr const char* NAME = "travel_time"; // The attribute's unique name.
};
