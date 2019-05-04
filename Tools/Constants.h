#pragma once

#include <limits>

// A special value representing infinity.
constexpr int INFTY = std::numeric_limits<int>::max() / 2;

// Special values representing an invalid (vertex/edge) ID.
constexpr int INVALID_ID = -1;
constexpr int INVALID_INDEX = -1;
constexpr int INVALID_VERTEX = -1;
constexpr int INVALID_EDGE = -1;

// This enum provides constants to specify the direction in which a road segment is open.
enum class RoadDirection { OPEN_IN_BOTH, FORWARD, REVERSE, CLOSED };

// The earth's mean radius in meters.
constexpr int EARTH_RADIUS = 6371000;

// The maximum number of shortest paths computed simultaneously during traffic assignment.
#ifndef TA_LOG_K
# define TA_LOG_K 5
#endif

// The maximum number of source vertices that are to be substituted into radiation model's formula.
#ifndef DC_MAX_NUM_SOURCES
# define DC_MAX_NUM_SOURCES std::numeric_limits<int>::max()
#endif
