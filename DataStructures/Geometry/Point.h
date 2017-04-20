#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <ostream>

// A point on a two-dimensional plane.
class Point {
 public:
  // Constructs a point.
  Point() : x(0), y(0) {}

  // Constructs a point with the specified coordinates.
  Point(const int x, const int y) : x(x), y(y) {}

  // Returns the x-coordinate.
  int getX() const {
    return x;
  }

  // Returns the y-coordinate.
  int getY() const {
    return y;
  }

  // Takes the coordinate-wise minimum of this and the specified point.
  void min(const Point& other) {
    x = std::min(x, other.x);
    y = std::min(y, other.y);
  }

  // Takes the coordinate-wise maximum of this and the specified point.
  void max(const Point& other) {
    x = std::max(x, other.x);
    y = std::max(y, other.y);
  }

  // Returns the Manhattan distance to the specified point.
  int getManhattanDistanceTo(const Point& other) const {
    return std::abs(x - other.x) + std::abs(y - other.y);
  }

  // Returns the Euclidean distance to the specified point.
  float getEuclideanDistanceTo(const Point& other) const {
    return std::sqrt(getSquaredEuclideanDistanceTo(other));
  }

  // Returns the squared Euclidean distance to the specified point.
  int64_t getSquaredEuclideanDistanceTo(const Point& other) const {
    int64_t deltaX = x - other.x;
    int64_t deltaY = y - other.y;
    return deltaX * deltaX + deltaY * deltaY;
  }

  // Some useful arithmetic operators.
  friend Point operator+(const Point& lhs, const Point& rhs) {
    return {lhs.getX() + rhs.getX(), lhs.getY() + rhs.getY()};
  }

  friend Point operator-(const Point& lhs, const Point& rhs) {
    return {lhs.getX() - rhs.getX(), lhs.getY() - rhs.getY()};
  }

  // Write a textual representation to the specified output stream.
  friend std::ostream& operator<<(std::ostream& os, const Point& point) {
    os << "(" << point.getX() << ", " << point.getY() << ")";
    return os;
  }

 private:
  int x; // The X coordinate.
  int y; // The Y coordinate.
};
