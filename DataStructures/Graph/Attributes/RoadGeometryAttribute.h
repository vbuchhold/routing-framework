#pragma once

#include <cassert>
#include <vector>

#include "DataStructures/Geometry/LatLng.h"
#include "DataStructures/Graph/Attributes/AbstractAttribute.h"

// An attribute associating road geometry data with each edge of a graph.
class RoadGeometryAttribute : public AbstractAttribute<std::vector<LatLng>> {
 public:
  // Returns a vector of LatLngs representing the road geometry of edge e.
  const Type& roadGeometry(const int e) const {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

  // Returns a reference to a vector of LatLngs representing the road geometry of edge e.
  Type& roadGeometry(const int e) {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

 protected:
  static constexpr const char* NAME = "road_geometry"; // The attribute's unique name.
};
