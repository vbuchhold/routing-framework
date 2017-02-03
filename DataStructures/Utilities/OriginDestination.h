#pragma once

// An origin-destination (O-D) pair, representing a travel demand or a query.
struct OriginDestination {
  // Constructs an O-D pair from origin to destination.
  OriginDestination(const int origin, const int destination)
      : origin(origin), destination(destination) {}

  int origin;
  int destination;
};
