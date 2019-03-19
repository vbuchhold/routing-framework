#pragma once

#include <cassert>

#include "DataStructures/Graph/Attributes/AbstractAttribute.h"

// An attribute associating a number of opportunities with each vertex of a graph.
class NumOpportunitiesAttribute : public AbstractAttribute<int> {
 public:
  // Returns the number of opportunities at vertex v.
  const Type& numOpportunities(const int v) const {
    assert(v >= 0); assert(v < values.size());
    return values[v];
  }

  // Returns a reference to the number of opportunities at vertex v.
  Type& numOpportunities(const int v) {
    assert(v >= 0); assert(v < values.size());
    return values[v];
  }

 protected:
  static constexpr const char* NAME = "num_opportunities"; // The attribute's unique name.
};
