#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <ostream>

#include "Tools/Constants.h"
#include "Tools/Math.h"
#include "Tools/Workarounds.h"

// This class represents a geographic point, specified by its latitude and longitude.
class LatLng {
 public:
  // Latitude and longitude are stored internally as integers in units of 1/PRECISION degrees.
  static constexpr int PRECISION = 1000000;
  static constexpr int DEG_90 = 90 * PRECISION;
  static constexpr int DEG_180 = 180 * PRECISION;
  static constexpr int DEG_360 = 360 * PRECISION;

  // Constructs an uninitialized LatLng.
  constexpr LatLng() : lat(INFTY), lng(INFTY) {}

  // Constructs a LatLng. Coordinates are specified in 1/PRECISION degrees. If wrap is set to true,
  // latitude is automatically clamped to the range [-90 deg, 90 deg] and longitude is
  // automatically wrapped so that it falls within the range [-180 deg, 180 deg].
  LatLng(const int lat, const int lng, const bool wrap = false) {
    // Clamp latitude.
    this->lat = wrap ? std::min(std::max(lat, -DEG_90), use(DEG_90)) : lat;

    // Wrap longitude.
    const int a = (lng + DEG_180) % DEG_360;
    this->lng = wrap ? (a < 0 ? a + DEG_180 : a - DEG_180) : lng;

	// Coordinates have to be valid at this point, even if automatic wrapping is disabled.
    assert(this->lat >= -DEG_90); assert(this->lat <= DEG_90);
    assert(this->lng >= -DEG_180); assert(this->lng <= DEG_180);
  }

  // Constructs a LatLng. Coordinates are specified in degrees. If wrap is set to true, latitude is
  // automatically clamped to the range [-90 deg, 90 deg] and longitude is automatically wrapped so
  // that it falls within the range [-180 deg, 180 deg].
  LatLng(const double lat, const double lng, const bool wrap = false)
      : LatLng(static_cast<int>(std::round(lat * PRECISION)),
               static_cast<int>(std::round(lng * PRECISION)), wrap) {}

  // Returns true if the LatLng represents a valid geographic point, false otherwise.
  bool isValid() const {
    return lat != INFTY;
  }

  // Returns the latitude in 1/PRECISION degrees.
  int latitude() const {
    assert(isValid());
    return lat;
  }

  // Returns the longitude in 1/PRECISION degrees.
  int longitude() const {
    assert(isValid());
    return lng;
  }

  // Returns the latitude in degrees.
  double latInDeg() const {
    assert(isValid());
    return lat / static_cast<double>(PRECISION);
  }

  // Returns the longitude in degrees.
  double lngInDeg() const {
    assert(isValid());
    return lng / static_cast<double>(PRECISION);
  }

  // Takes the coordinate-wise minimum of this and the specified LatLng.
  void min(const LatLng& other) {
    assert(isValid()); assert(other.isValid());
    lat = std::min(lat, other.lat);
    lng = std::min(lng, other.lng);
  }

  // Takes the coordinate-wise maximum of this and the specified LatLng.
  void max(const LatLng& other) {
    assert(isValid()); assert(other.isValid());
    lat = std::max(lat, other.lat);
    lng = std::max(lng, other.lng);
  }

  // Returns the great-circle distance in meters to the specified LatLng using haversine formula.
  // For more details, see: www.movable-type.co.uk/scripts/latlong.html
  double getGreatCircleDistanceTo(const LatLng& other) const {
    assert(isValid()); assert(other.isValid());
    const double lat1 = toRadians(latInDeg());
    const double lat2 = toRadians(other.latInDeg());
    const double deltaLat = toRadians(other.latInDeg() - latInDeg());
    const double deltaLng = toRadians(other.lngInDeg() - lngInDeg());

    const double a = std::sin(deltaLat / 2) * std::sin(deltaLat / 2) +
        std::cos(lat1) * std::cos(lat2) * std::sin(deltaLng / 2) * std::sin(deltaLng / 2);
    const double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));

    return EARTH_RADIUS * c;
  }

 private:
  int lat; // The latitude in 1/PRECISION degrees.
  int lng; // The longitude in 1/PRECISION degrees.
};

// Some useful arithmetic operators. Automatic wrapping takes place.
inline LatLng operator+(const LatLng& lhs, const LatLng& rhs) {
  assert(lhs.isValid()); assert(rhs.isValid());
  return LatLng(lhs.latitude() + rhs.latitude(), lhs.longitude() + rhs.longitude(), true);
}

inline LatLng operator-(const LatLng& lhs, const LatLng& rhs) {
  assert(lhs.isValid()); assert(rhs.isValid());
  return LatLng(lhs.latitude() - rhs.latitude(), lhs.longitude() - rhs.longitude(), true);
}

// Write a textual representation to the specified output stream.
inline std::ostream& operator<<(std::ostream& os, const LatLng& latLng) {
  os << "(" << latLng.latitude() << ", " << latLng.longitude() << ")";
  return os;
}
