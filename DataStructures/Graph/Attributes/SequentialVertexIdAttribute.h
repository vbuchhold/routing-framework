#pragma once

#include <cassert>

#include "DataStructures/Graph/Attributes/AbstractAttribute.h"
#include "Tools/Constants.h"

// An attribute associating a sequential original ID with each vertex of a graph. The attribute can
// be used to store, for each vertex of a subgraph, the corresponding vertex in the supergraph.
class SequentialVertexIdAttribute : public AbstractAttribute<int> {
 public:
  // Returns the attribute's default value.
  static Type defaultValue() {
    return INVALID_VERTEX;
  }

  // Returns a reference to the sequential original ID of vertex v.
  Type& sequentialVertexId(const int v) {
    assert(v >= 0); assert(v < values.size());
    return values[v];
  }

  // Returns a const reference to the sequential original ID of vertex v.
  const Type& sequentialVertexId(const int v) const {
    assert(v >= 0); assert(v < values.size());
    return values[v];
  }

 protected:
  static constexpr const char* NAME = "sequential_vertex_id"; // The attribute's unique name.
};
