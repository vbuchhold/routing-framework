#pragma once

#include <cassert>

#include "DataStructures/Graph/Attributes/AbstractAttribute.h"
#include "Tools/Constants.h"

// An attribute associating an ID with each vertex of a graph.
class VertexIdAttribute : public AbstractAttribute<int> {
 public:
  static constexpr Type DEFAULT_VALUE = INVALID_ID; // The attribute's default value.
  static constexpr const char* NAME = "vertex_id";  // The attribute's unique name.

  // Returns the ID of vertex v.
  const Type& vertexId(const int v) const {
    assert(v >= 0); assert(v < values.size());
    return values[v];
  }

  // Returns a reference to the ID of vertex v.
  Type& vertexId(const int v) {
    assert(v >= 0); assert(v < values.size());
    return values[v];
  }
};
