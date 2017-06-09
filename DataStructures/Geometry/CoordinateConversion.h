#pragma once

#include <cassert>
#include <cmath>
#include <string>

#include <proj_api.h>

#include "DataStructures/Geometry/LatLng.h"
#include "DataStructures/Geometry/Point.h"
#include "Tools/Math.h"

// This class provides coordinate conversion between a primary coordinate system and the World
// Geodetic System 1984 (WGS84). Since the class uses the proj.4 library, the primary coordinate
// system has to be specified as a proj.4 string (or as an EPSG code).
// For more details, see: www.proj4.org
class CoordinateConversion {
 public:
  // proj.4 strings for some widely used coordinate systems.
  static constexpr const char* WGS84 = // EPSG::4326
      "+proj=longlat +datum=WGS84 +no_defs";
  static constexpr const char* DHDN_GAUSS_KRUGER_ZONE_3 = // EPSG::31467
      "+proj=tmerc +lat_0=0 +lon_0=9 +k=1 +x_0=3500000 +y_0=0 +datum=potsdam +no_defs";

  // Constructs a conversion between WGS84 and the coordinate system specified as a proj.4 string.
  CoordinateConversion(const std::string& primaryCrs)
      : wgs84Crs(pj_init_plus(WGS84)), otherCrs(pj_init_plus(primaryCrs.c_str())) {
    assert(wgs84Crs);
    assert(otherCrs);
  }

  // Constructs a conversion between WGS84 and the coordinate system specified as an EPSG code.
  CoordinateConversion(const int primaryCrs)
      : wgs84Crs(pj_init_plus(WGS84)),
        otherCrs(pj_init_plus(("+init=epsg:" + std::to_string(primaryCrs)).c_str())) {
    assert(wgs84Crs);
    assert(otherCrs);
  }

  // Releases all resources.
  ~CoordinateConversion() {
    pj_free(wgs84Crs);
    pj_free(otherCrs);
  }

  // Returns the coordinate given in the primary coordinate system as a WGS84 coordinate.
  LatLng convert(const Point& p) {
    double x = p.getX();
    double y = p.getY();
    if (pj_transform(otherCrs, wgs84Crs, 1, 1, &x, &y, nullptr))
      assert(false);
    assert(x != HUGE_VAL);
    assert(y != HUGE_VAL);
    return LatLng(toDegrees(y), toDegrees(x));
  }

  // Returns the WGS84 coordinate as a coordinate in the primary coordinate system.
  Point convert(const LatLng& l) {
    double x = toRadians(l.lngInDeg());
    double y = toRadians(l.latInDeg());
    if (pj_transform(wgs84Crs, otherCrs, 1, 1, &x, &y, nullptr))
      assert(false);
    assert(x != HUGE_VAL);
    assert(y != HUGE_VAL);
    return Point(std::round(x), std::round(y));
  }

 private:
  const projPJ wgs84Crs; // The WGS84 coordinate system (CRS).
  const projPJ otherCrs; // The primary coordinate system (CRS).
};
