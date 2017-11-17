#pragma once

#include <algorithm>
#include <cassert>
#include <ostream>
#include <vector>

#include "DataStructures/Geometry/Helpers.h"
#include "DataStructures/Geometry/Point.h"

// A polygon defines a two-dimensional region enclosed by a single closed polygonal chain.
class Polygon {
 public:
  // Iterators referring to vertices of this polygon.
  using VertexIterator = std::vector<Point>::iterator;
  using ConstVertexIterator = std::vector<Point>::const_iterator;
  using ReverseVertexIterator = std::vector<Point>::reverse_iterator;
  using ConstReverseVertexIterator = std::vector<Point>::const_reverse_iterator;

  // Constructs an empty polygon.
  Polygon() = default;

  // Constructs a polygon with vertices from the specified sequence.
  template <typename VertexIteratorT>
  Polygon(VertexIteratorT first, VertexIteratorT last) {
    add(first, last);
  }

  // Writes a character representation to the specified output stream.
  friend std::ostream& operator<<(std::ostream& os, const Polygon& polygon) {
    os << "(";
    for (int i = 0; i < polygon.size(); ++i) {
      if (i != 0)
        os << ", ";
      os << polygon.vertices[i];
    }
    os << ")";
    return os;
  }

  // Returns the vertex with the specified index.
  const Point& operator[](const int idx) const {
    assert(idx >= 0); assert(idx < size());
    return vertices[idx];
  }

  // Returns true if this polygon is of size 0.
  bool empty() const {
    return vertices.empty();
  }

  // Returns the number of vertices.
  int size() const {
    return vertices.size();
  }

  // Returns the vertex with the lowest index.
  const Point& front() const {
    assert(!empty());
    return vertices.front();
  }

  // Returns the vertex with the highest index.
  const Point& back() const {
    assert(!empty());
    return vertices.back();
  }

  // Returns an iterator referring to the first vertex of this polygon.
  ConstVertexIterator begin() const {
    return vertices.begin();
  }

  // Returns an iterator which is the past-the-end value for this polygon.
  ConstVertexIterator end() const {
    return vertices.end();
  }

  ConstReverseVertexIterator rbegin() const {
    return vertices.rbegin();
  }

  ConstReverseVertexIterator rend() const {
    return vertices.rend();
  }

  // Appends the specified vertex to this polygon.
  void add(const Point& p) {
    vertices.push_back(p);
  }

  // Appends all of the specified vertices to this polygon.
  template <typename VertexIteratorT>
  void add(VertexIteratorT first, VertexIteratorT last) {
    vertices.insert(vertices.end(), first, last);
  }

  // Removes the vertex with the highest index from this polygon.
  void removeBack() {
    assert(!empty());
    vertices.pop_back();
  }

  // Removes the vertex pointed to by pos from this polygon.
  VertexIterator remove(ConstVertexIterator pos) {
    assert(vertices.begin() <= pos); assert(pos < vertices.end());
    return vertices.erase(pos);
  }

  // Returns the index of a lowest of the leftmost vertices.
  // Precondition: The polygon must not be empty.
  int leftmostVertex() const {
    assert(!empty());
    int minIdx = 0;
    Point min = vertices[0];
    for (int i = 1; i < size(); ++i)
      if (vertices[i].getX() < min.getX() ||
          (vertices[i].getX() == min.getX() && vertices[i].getY() < min.getY())) {
        minIdx = i;
        min = vertices[i];
      }
    return minIdx;
  }

  // Returns 1 or -1 as this polygon is counterclockwise or clockwise oriented.
  // Precondition: The polygon must be simple.
  int orientation() const {
    assert(simple());
    const int q = leftmostVertex();
    const int p = (q - 1 + size()) % size();
    const int r = (q + 1) % size();
    return ::orientation(vertices[p], vertices[q], vertices[r]);
  }

  // Reverses the orientation of this polygon.
  void reverseOrientation() {
    std::reverse(vertices.begin(), vertices.end());
  }

  // Returns true if this polygon is simple, i.e., its boundary does not intersect itself.
  // Note: The implementation takes quadratic time, whereas the problem can be solved in O(nlog⁡n).
  bool simple() const {
    if (size() < 3)
      return false;
    // Test each pair of edges. Edges that are consecutive must not form an angle of zero degrees.
    // Non-consecutive edges must not intersect.
    for (int i = size() - 1, j = 0; j < size() - 1; i = j++) {
      // The consecutive edges ij and j(j + 1) must not form an angle of zero degrees.
      if (::orientation(vertices[i], vertices[j], vertices[j + 1]) == 0 &&
          (vertices[i] - vertices[j]) * (vertices[j + 1] - vertices[j]) > 0)
        return false;
      for (int k = j + 1; k < size() - 1; ++k)
        if (k + 1 != i) {
          // The non-consecutive edges ij and k(k + 1) must not intersect.
          if (intersection(vertices[i], vertices[j], vertices[k], vertices[k + 1]))
            return false;
        } else {
          // The consecutive edges ki and ij must not form an angle of zero degrees.
          if (::orientation(vertices[k], vertices[i], vertices[j]) == 0 &&
              (vertices[k] - vertices[i]) * (vertices[j] - vertices[i]) > 0)
            return false;
        }
    }
    return true;
  }

  // Returns true if r is inside the boundary of this polygon.
  bool contains(const Point& r) const {
    bool inside = false;
    for (int i = size() - 1, j = 0; j < size(); i = j++)
      if ((r.getY() < vertices[i].getY()) != (r.getY() < vertices[j].getY()) &&
          ::orientation(vertices[i], vertices[j], r) * (vertices[j].getY()-vertices[i].getY()) > 0)
        inside = !inside;
    return inside;
  }

 private:
  std::vector<Point> vertices; // The vertices of this polygon.
};
