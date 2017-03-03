#pragma once

#include <cassert>
#include <cstdint>

#include <boost/dynamic_bitset.hpp>

#include "DataStructures/Geometry/PointSet.h"

namespace kdtree {

// An adapter that makes a subset of a PointSet searchable by a kd-tree.
class PointSubsetAdapter {
 public:
  // Constructs an adapter for the specified subset of the given PointSet.
  PointSubsetAdapter(const PointSet& pointSet,
                     const int subsetSize, const boost::dynamic_bitset<>& subset)
      : pointSet(pointSet) {
    updateSubset(subsetSize, subset);
  }

  // Returns the number of points in the whole PointSet (not only in the subset).
  int numPointsInPointSet() const {
    return pointSet.size();
  }

  // Returns the point with index idx.
  const Point& operator[](const int idx) const {
    return pointSet[idx];
  }

  // Returns the index in the PointSet of the point with index idx.
  int getIndexInPointSet(const int idx) const {
    assert(idx >= 0); assert(idx < indexInPointSet.size());
    return indexInPointSet[idx];
  }

  // Updates the subset that should be searched.
  void updateSubset(const int subsetSize, const boost::dynamic_bitset<>& subset) {
    assert(subsetSize >= 0);
    indexInPointSet.resize(subsetSize);
    int nextIdxInSubset = 0;
    for (int i = subset.find_first(); i < boost::dynamic_bitset<>::npos; i = subset.find_next(i)) {
      assert(nextIdxInSubset < indexInPointSet.size());
      indexInPointSet[nextIdxInSubset++] = i;
    }
  }

  // Returns the number of points in the subset.
  int kdtree_get_point_count() const {
    return indexInPointSet.size();
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

  // Tell a kd-tree to compute the bounding box.
  template <typename BoundingBoxT>
  bool kdtree_get_bbox(BoundingBoxT& /*bb*/) const {
    return false;
  }

 private:
  const PointSet& pointSet;         // The set of points of which a subset should be searched.
  std::vector<int> indexInPointSet; // For each point p in the subset, p's index in the PointSet.
};

}
