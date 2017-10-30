#pragma once

#include <cassert>

#include "DataStructures/Geometry/Point.h"
#include "DataStructures/Graph/Attributes/AbstractAttribute.h"

// An attribute associating a coordinate with each vertex of a graph.
class CoordinateAttribute : public AbstractAttribute<Point> {
 public:
  // Returns the coordinate of vertex v.
  const Type& coordinate(const int v) const {
    assert(v >= 0); assert(v < values.size());
    return values[v];
  }

  // Returns a reference to the coordinate of vertex v.
  Type& coordinate(const int v) {
    assert(v >= 0); assert(v < values.size());
    return values[v];
  }

 protected:
  static constexpr const char* NAME = "coordinate"; // The attribute's unique name.
};
