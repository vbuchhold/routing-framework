#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <ostream>

// A point on a two-dimensional plane.
class Point {
 public:
  // Constructs a point at the origin.
  Point() : x(0), y(0) {}

  // Constructs a point at the specified location.
  Point(const int x, const int y) : x(x), y(y) {}

  // Some useful arithmetic operators.
  friend Point operator+(const Point& lhs, const Point& rhs) {
    return {lhs.x + rhs.x, lhs.y + rhs.y};
  }

  friend Point operator-(const Point& lhs, const Point& rhs) {
    return {lhs.x - rhs.x, lhs.y - rhs.y};
  }

  friend int64_t operator*(const Point& lhs, const Point& rhs) {
    return static_cast<int64_t>(lhs.x) * rhs.x + static_cast<int64_t>(lhs.y) * rhs.y;
  }

  // Writes a character representation to the specified output stream.
  friend std::ostream& operator<<(std::ostream& os, const Point& p) {
    os << "(" << p.x << ", " << p.y << ")";
    return os;
  }

  // Returns true if lhs and rhs coincide.
  friend bool operator==(const Point& lhs, const Point& rhs) {
    return lhs.x == rhs.x && lhs.y == rhs.y;
  }

  // Returns true if lhs and rhs do not coincide.
  friend bool operator!=(const Point& lhs, const Point& rhs) {
    return !(lhs == rhs);
  }

  // Returns the x-coordinate.
  int getX() const {
    return x;
  }

  // Returns the y-coordinate.
  int getY() const {
    return y;
  }

  // Returns a reference to the x-coordinate.
  int& getX() {
    return x;
  }

  // Returns a reference to the y-coordinate.
  int& getY() {
    return y;
  }

  // Takes the coordinate-wise minimum of this and the specified point.
  void min(const Point& p) {
    x = std::min(x, p.x);
    y = std::min(y, p.y);
  }

  // Takes the coordinate-wise maximum of this and the specified point.
  void max(const Point& p) {
    x = std::max(x, p.x);
    y = std::max(y, p.y);
  }

  // Returns the Manhattan or "city block" distance to the specified point.
  int getManhattanDistanceTo(const Point& p) const {
    return std::abs(x - p.x) + std::abs(y - p.y);
  }

  // Returns the Euclidean distance to the specified point.
  double getEuclideanDistanceTo(const Point& p) const {
    return std::sqrt(getSquaredEuclideanDistanceTo(p));
  }

  // Returns the squared Euclidean distance to the specified point.
  int64_t getSquaredEuclideanDistanceTo(const Point& p) const {
    int64_t deltaX = x - p.x;
    int64_t deltaY = y - p.y;
    return deltaX * deltaX + deltaY * deltaY;
  }

  // Returns the chessboard, Chebyshev, or sup norm distance to the specified point.
  int getChebyshevDistanceTo(const Point& p) const {
    return std::max(std::abs(x - p.x), std::abs(y - p.y));
  }

 private:
  int x; // The X coordinate.
  int y; // The Y coordinate.
};
