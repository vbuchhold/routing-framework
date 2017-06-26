#pragma once

#include <cassert>
#include <ostream>

#include "DataStructures/Geometry/Point.h"

// A rectangle on a two-dimensional plane with sides parallel to the x- and y-axis.
class Rectangle {
 public:
  // Constructs a rectangle from a single point.
  explicit Rectangle(const Point& p = Point()) : southWest(p), northEast(p) {}

  // Constructs a rectangle from the points at its south-west and north-east corners.
  Rectangle(const Point& southWest, const Point& northEast)
      : southWest(southWest), northEast(northEast) {
    assert(southWest.getX() <= northEast.getX());
    assert(southWest.getY() <= northEast.getY());
  }

  // Constructs a bounding box containing all specified points.
  template <typename PointIteratorT>
  Rectangle(PointIteratorT first, PointIteratorT last) {
    if (first != last) {
      southWest = *first;
      northEast = *first;
    }
    extend(++first, last);
  }

  // Writes a character representation to the specified output stream.
  friend std::ostream& operator<<(std::ostream& os, const Rectangle& rect) {
    os << "(SW=" << rect.southWest << ", NE=" << rect.northEast << ")";
    return os;
  }

  // Returns the south-west corner.
  const Point& getSouthWest() const {
    return southWest;
  }

  // Returns the north-east corner.
  const Point& getNorthEast() const {
    return northEast;
  }

  // Returns true if p is inside the boundary of this rectangle.
  bool contains(const Point& p) const {
    return southWest.getX() <= p.getX() && p.getX() <= northEast.getX() &&
        southWest.getY() <= p.getY() && p.getY() <= northEast.getY();
  }

  // Returns true if this rectangle and the specified rectangle intersect.
  bool intersects(const Rectangle& rect) const {
    return northEast.getX() >= rect.southWest.getX() && southWest.getX() <= rect.northEast.getX() &&
        northEast.getY() >= rect.southWest.getY() && southWest.getY() <= rect.northEast.getY();
  }

  // Extends this rectangle to contain the point p.
  void extend(const Point& p) {
    southWest.min(p);
    northEast.max(p);
  }

  // Extends this rectangle to contain the specified points.
  template <typename PointIteratorT>
  void extend(PointIteratorT first, PointIteratorT last) {
    for (; first != last; ++first)
      extend(*first);
  }

 private:
  Point southWest; // The south-west corner.
  Point northEast; // The north-east corner.
};
