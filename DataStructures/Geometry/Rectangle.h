#pragma once

#include <cassert>

#include "DataStructures/Geometry/Point.h"

// A rectangle on a two-dimensional plane.
class Rectangle {
 public:
  // Constructs a rectangle from a single point.
  explicit Rectangle(const Point& p = Point()) : southWest(p), northEast(p) {}

  // Constructs a rectangle from two points.
  Rectangle(const Point& southWest, const Point& northEast)
      : southWest(southWest), northEast(northEast) {
    assert(southWest.getX() <= northEast.getX());
    assert(southWest.getY() <= northEast.getY());
  }

  // Returns the south-west corner.
  const Point& getSouthWest() const {
    return southWest;
  }

  // Returns the north-east corner.
  const Point& getNorthEast() const {
    return northEast;
  }

  // Return true if this rectangle encloses the point p.
  bool contains(const Point& p) const {
    return southWest.getX() <= p.getX() && p.getX() <= northEast.getX() &&
        southWest.getY() <= p.getY() && p.getY() <= northEast.getY();
  }

  // Extends this rectangle to contain the point p.
  void extend(const Point& p) {
    southWest.min(p);
    northEast.max(p);
  }

 private:
  Point southWest; // The south-west corner.
  Point northEast; // The north-east corner.
};
