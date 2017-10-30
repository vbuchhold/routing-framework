#pragma once

#include <cassert>

#include "DataStructures/Graph/Attributes/AbstractAttribute.h"

// Road categories defined by the XATF file format.
enum class XatfRoadCategory {
  MOTORWAY_FAST        = 1,
  MOTORWAY_MEDIUM      = 2,
  MOTORWAY_SLOW        = 3,
  NATIONAL_ROAD_FAST   = 4,
  NATIONAL_ROAD_MEDIUM = 5,
  NATIONAL_ROAD_SLOW   = 6,
  REGIONAL_ROAD_FAST   = 7,
  REGIONAL_ROAD_MEDIUM = 8,
  REGIONAL_ROAD_SLOW   = 9,
  URBAN_STREET_FAST    = 10,
  URBAN_STREET_MEDIUM  = 11,
  URBAN_STREET_SLOW    = 12,
  FERRY                = 13,
  UNUSED               = 14,
  FOREST_ROAD          = 15,
};

// An attribute associating an XATF road category with each edge of a graph.
class XatfRoadCategoryAttribute : public AbstractAttribute<XatfRoadCategory> {
 public:
  // Returns the attribute's default value.
  static Type defaultValue() {
    return XatfRoadCategory::UNUSED;
  }

  // Returns the XATF road category of edge e.
  const Type& xatfRoadCategory(const int e) const {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

  // Returns a reference to the XATF road category of edge e.
  Type& xatfRoadCategory(const int e) {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

 protected:
  static constexpr const char* NAME = "xatf_road_category"; // The attribute's unique name.
};
