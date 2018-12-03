#pragma once

#include <cassert>

#include "DataStructures/Graph/Attributes/AbstractAttribute.h"

// An attribute associating a population size with each vertex of a graph.
class PopulationAttribute : public AbstractAttribute<int> {
 public:
  // Returns the population of vertex v.
  const Type& population(const int v) const {
    assert(v >= 0); assert(v < values.size());
    return values[v];
  }

  // Returns a reference to the population of vertex v.
  Type& population(const int v) {
    assert(v >= 0); assert(v < values.size());
    return values[v];
  }

 protected:
  static constexpr const char* NAME = "population"; // The attribute's unique name.
};
