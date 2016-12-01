#pragma once

#include <cassert>

#include "DataStructures/Geometry/LatLng.h"
#include "DataStructures/Graph/Attributes/AbstractAttribute.h"

// An attribute associating a LatLng with each vertex of a graph.
class LatLngAttribute : public AbstractAttribute<LatLng> {
 public:
  // A functor that returns the LatLng of the specified vertex in the specified graph. Used for
  // telling algorithms on which attribute of a graph they should work.
  template <typename GraphT>
  struct GetLatLng {
    const Type& operator()(const GraphT& graph, const int e) const { return graph.latLng(e); }
    Type& operator()(GraphT& graph, const int e) const { return graph.latLng(e); }
  };

  static constexpr const Type& DEFAULT_VALUE = Type();    // The attribute's default value.
  static constexpr const char* NAME          = "lat_lng"; // The attribute's unique name.

  // Returns a const reference to the LatLng of vertex v.
  const Type& latLng(const int v) const {
    assert(v >= 0); assert(v < values.size());
    return values[v];
  }

  // Returns a reference to the LatLng of vertex v.
  Type& latLng(const int v) {
    assert(v >= 0); assert(v < values.size());
    return values[v];
  }
};
