#pragma once

#include <cassert>

#include "DataStructures/Graph/Attributes/AbstractAttribute.h"
#include "Tools/Constants.h"

// An attribute associating a generic traversal cost with each edge of a graph.
class TraversalCostAttribute : public AbstractAttribute<int> {
 public:
  // Returns the attribute's default value.
  static Type defaultValue() {
    return INFTY;
  }

  // Returns the traversal cost of edge e.
  const Type& traversalCost(const int e) const {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

  // Returns a reference to the traversal cost of edge e.
  Type& traversalCost(const int e) {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

 protected:
  static constexpr const char* NAME = "traversal_cost"; // The attribute's unique name.
};
