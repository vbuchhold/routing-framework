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

// Sets min to the smaller and max to the greater of a and b.
template <typename T>
inline constexpr void minmax(const T a, const T b, T& min, T& max) {
  if (b < a) {
    min = b;
    max = a;
  } else {
    min = a;
    max = b;
  }
}

// Returns 0, 1, or -1 as the specified value is equal to, greater than, or less than zero.
template <typename T>
inline int signum(const T val) {
  return (0 < val) - (val < 0);
}
