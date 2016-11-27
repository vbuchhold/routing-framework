#pragma once

#include <cassert>

#include "DataStructures/Geometry/LatLng.h"
#include "DataStructures/Graph/Attributes/AbstractAttribute.h"

// An attribute associating a LatLng with each vertex of a graph.
class LatLngAttribute : public AbstractAttribute<LatLng> {
 public:
  static constexpr const Type& DEFAULT_VALUE = Type();    // The attribute's default value.
  static constexpr const char* NAME          = "lat_lng"; // The attribute's unique name.

  // Returns the LatLng of vertex v.
  Type latLng(const int v) const {
    assert(v >= 0); assert(v < values.size());
    return values[v];
  }

  // Sets the LatLng of vertex v to the specified value.
  void setLatLng(const int v, const Type& val) {
    assert(v >= 0); assert(v < values.size());
    assert(val.isValid());
    values[v] = val;
  }
};
