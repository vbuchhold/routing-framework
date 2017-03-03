#pragma once

#include <cassert>
#include <vector>

#include "DataStructures/Geometry/Rectangle.h"

// A set of points on a two-dimensional plane.
class PointSet {
 public:
  // Constructs an empty set of points.
  explicit PointSet(const int minCapacity = 0) {
    reserve(minCapacity);
  }

  // Returns the point with index idx.
  const Point& operator[](const int idx) const {
    assert(idx >= 0); assert(idx < points.size());
    return points[idx];
  }

  // Returns a bounding box containing all points of this set.
  const Rectangle& getBoundingBox() const {
    assert(!empty());
    return boundingBox;
  }

  // Returns true if this set contains no points.
  bool empty() const {
    return points.empty();
  }

  // Returns the number of points in this set.
  int size() const {
    return points.size();
  }

  // Ensures that this set can hold at least minCapacity points without requiring reallocation.
  void reserve(const int minCapacity) {
    assert(minCapacity >= 0);
    points.reserve(minCapacity);
  }

  // Inserts the point p into this set.
  void insert(const Point& p) {
    if (points.empty())
      boundingBox = Rectangle(p);
    else
      boundingBox.extend(p);
    points.push_back(p);
  }

  // Removes all points from this set.
  void clear() {
    points.clear();
  }

 private:
  std::vector<Point> points; // The points of this set.
  Rectangle boundingBox;     // A bounding box containing all points of this set.
};
