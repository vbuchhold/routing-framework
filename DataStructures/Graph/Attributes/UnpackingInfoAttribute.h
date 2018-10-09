#pragma once

#include <cassert>

#include "DataStructures/Graph/Attributes/AbstractAttribute.h"
#include "Tools/Constants.h"

// An attribute associating unpacking information with each edge of a CH search graph. The unpacking
// information is a pair of 32-bit integers. In the case of shortcut edges, the first integer is the
// index of the first originating edge in the downward graph, and the second integer is the index of
// the second originating edge in the upward graph. In the case of input edges, the first integer is
// the index of this edge in the input graph, and the second integer is set to INVALID_EDGE.
class UnpackingInfoAttribute : public AbstractAttribute<std::pair<int, int>> {
 public:
  // Returns the attribute's default value.
  static Type defaultValue() {
    return {INVALID_EDGE, INVALID_EDGE};
  }

  // Returns the unpacking information of edge e.
  const Type& unpackingInfo(const int e) const {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

  // Returns a reference to the unpacking information of edge e.
  Type& unpackingInfo(const int e) {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

 protected:
  static constexpr const char* NAME = "unpacking_info"; // The attribute's unique name.
};
