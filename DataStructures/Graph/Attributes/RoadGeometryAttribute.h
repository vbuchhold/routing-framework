#pragma once

#include <cassert>
#include <vector>

#include "DataStructures/Geometry/LatLng.h"
#include "DataStructures/Graph/Attributes/AbstractAttribute.h"

// An attribute associating road geometry data with each edge of a graph.
class RoadGeometryAttribute : public AbstractAttribute<std::vector<LatLng>> {
 public:
  // Returns a vector of LatLngs representing the road geometry of edge e.
  const Type& roadGeometry(const int v) const {
    assert(v >= 0); assert(v < values.size());
    return values[v];
  }

  // Returns a reference to a vector of LatLngs representing the road geometry of edge e.
  Type& roadGeometry(const int v) {
    assert(v >= 0); assert(v < values.size());
    return values[v];
  }

 protected:
  static constexpr const char* NAME = "road_geometry"; // The attribute's unique name.
};
