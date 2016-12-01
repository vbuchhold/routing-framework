#pragma once

#include <limits>

// A special value representing infinity.
constexpr int INFTY = std::numeric_limits<int>::max() / 2;

// Special values representing an invalid vertex/edge ID.
constexpr int INVALID_VERTEX = -1;
constexpr int INVALID_EDGE = -1;

// The earth's mean radius in meters.
constexpr int EARTH_RADIUS = 6371000;
