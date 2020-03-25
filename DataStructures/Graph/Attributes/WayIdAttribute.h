#pragma once

#include <cassert>

#include "DataStructures/Graph/Attributes/AbstractAttribute.h"
#include "Tools/Constants.h"

// An attribute storing explicitly the way id for each edge of graph
class WayIdAttribute : public AbstractAttribute<int> {
 public:
  // Returns the way id of edge e.
  const Type& wayId(const int e) const {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

  // Returns a reference to the way id of edge e.
  Type& wayId(const int e) {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

 protected:
  static constexpr const char* NAME = "way_id"; // The attribute's unique name.
};
