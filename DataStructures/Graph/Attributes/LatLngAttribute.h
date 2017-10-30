#pragma once

#include <cassert>

#include "DataStructures/Geometry/LatLng.h"
#include "DataStructures/Graph/Attributes/AbstractAttribute.h"

// An attribute associating a LatLng with each vertex of a graph.
class LatLngAttribute : public AbstractAttribute<LatLng> {
 public:
  // Returns the LatLng of vertex v.
  const Type& latLng(const int v) const {
    assert(v >= 0); assert(v < values.size());
    return values[v];
  }

  // Returns a reference to the LatLng of vertex v.
  Type& latLng(const int v) {
    assert(v >= 0); assert(v < values.size());
    return values[v];
  }

 protected:
  static constexpr const char* NAME = "lat_lng"; // The attribute's unique name.
};
