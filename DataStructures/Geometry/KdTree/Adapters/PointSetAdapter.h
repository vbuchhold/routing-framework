#pragma once

#include <cassert>
#include <cstdint>

#include "DataStructures/Geometry/PointSet.h"
#include "DataStructures/Geometry/Rectangle.h"

namespace kdtree {

// An adapter that makes a PointSet searchable using a kd-tree.
class PointSetAdapter {
 public:
  // Constructs an adapter for the specified PointSet.
  PointSetAdapter(const PointSet& pointSet) : pointSet(pointSet) {}

  // Returns the point with index idx.
  const Point& operator[](const int idx) const {
    return pointSet[idx];
  }

  // Returns the index in the PointSet of the point with index idx.
  int getIndexInPointSet(const int idx) const {
    assert(idx >= 0); assert(idx < pointSet.size());
    return idx;
  }

  // Returns the number of points in the set.
  int kdtree_get_point_count() const {
    return pointSet.size();
  }

  // Returns the specified coordinate of the point with index idx.
  int kdtree_get_pt(const int idx, const int coordinate) const {
    assert(coordinate >= 0); assert(coordinate < 2);
    const int idxInPointSet = getIndexInPointSet(idx);
    return coordinate == 0 ? pointSet[idxInPointSet].getX() : pointSet[idxInPointSet].getY();
  }

  // Returns the squared Euclidean distance between the point p1 and the point with index p2Idx.
  int64_t kdtree_distance(const int* const p1, const int p2Idx, const int /*dimension*/) const {
    return pointSet[getIndexInPointSet(p2Idx)].getSquaredEuclideanDistanceTo({p1[0], p1[1]});
  }

  // Sets the bounding box bb of a kd-tree.
  template <typename BoundingBoxT>
  bool kdtree_get_bbox(BoundingBoxT& bb) const {
    const Rectangle& boundingBox = pointSet.getBoundingBox();
    bb[0].low = boundingBox.getSouthWest().getX();
    bb[1].low = boundingBox.getSouthWest().getY();
    bb[0].high = boundingBox.getNorthEast().getX();
    bb[1].high = boundingBox.getNorthEast().getY();
    return true;
  }

 private:
  const PointSet& pointSet; // The set of points that should be searched.
};

}
