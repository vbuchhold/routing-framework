#pragma once

// The ratio of the circumference of a circle to its diameter.
constexpr double PI = 3.14159265358979323846l;

// Converts an angle measured in degrees to an approximately equivalent angle measured in radians.
inline constexpr double toRadians(const double angdeg) {
  return PI / 180 * angdeg;
}

// Converts an angle measured in radians to an approximately equivalent angle measured in degrees.
inline constexpr double toDegrees(const double angrad) {
  return 180 / PI * angrad;
}
