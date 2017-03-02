#pragma once

#include <cassert>

#include <boost/dynamic_bitset.hpp>

#include "DataStructures/Geometry/KdTree/Adapters/PointSubsetAdapter.h"
#include "DataStructures/Geometry/KdTree/Metrics/EuclideanMetric.h"
#include "DataStructures/Geometry/KdTree/StaticKdTree.h"

// A dynamic kd-tree for sets of two-dimensional points. It supports efficient circular range
// queries and k-nearest-neighbor finding. The metric used to compute distances (Euclidean or
// Manhattan) is specified by a template parameter.
template <
    typename PointSetAdapterT = kdtree::PointSubsetAdapter,
    template <typename> class MetricT = kdtree::EuclideanMetric>
class DynamicKdTree {
 private:
  using Distance = typename MetricT<PointSetAdapterT>::DistanceType; // The distance type.

 public:
  // Constructs a dynamic kd-tree using the specified PointSet adapter.
  DynamicKdTree(PointSetAdapterT& pointSet) : pointSet(pointSet), kdTree(pointSet) {
    rebuild();
  }

  // Rebuilds this kd-tree.
  void rebuild() {
    kdTree.rebuild();
    isValidPoint.clear();
    isValidPoint.resize(pointSet.numPointsInPointSet());
    for (int i = 0; i < pointSet.kdtree_get_point_count(); ++i)
      isValidPoint[pointSet.getIndexInPointSet(i)] = true;
    numValidPoints = pointSet.kdtree_get_point_count();
  }

  // Returns true if this kd-tree contains the point with index idx.
  bool contains(const int idx) const {
    assert(idx >= 0); assert(idx < isValidPoint.size());
    return isValidPoint[idx];
  }

  // Invokes for each point within the specified circular query range the given hook function.
  template <typename PointFoundT>
  void circularRangeQuery(const Point& p, const Distance radius, PointFoundT pointFound) const {
    kdTree.circularRangeQuery(p, radius, [this, pointFound](const int idx, const Distance dist) {
      assert(idx >= 0); assert(idx < isValidPoint.size());
      if (isValidPoint[idx])
        pointFound(idx, dist);
    });
  }

  // Removes the point with index idx from this kd-tree.
  void remove(const int idx) {
    assert(idx >= 0); assert(idx < isValidPoint.size());
    assert(isValidPoint[idx]);
    isValidPoint[idx] = false;
    --numValidPoints;

    // Whenever the number of valid points becomes half the size of the kd-tree, we rebuild it.
    if (2 * numValidPoints <= pointSet.kdtree_get_point_count()) {
      pointSet.updateSubset(numValidPoints, isValidPoint);
      kdTree.rebuild();
    }
  }

 private:
  PointSetAdapterT& pointSet;                     // An adapter to the set of points.
  StaticKdTree<PointSetAdapterT, MetricT> kdTree; // The underlying static kd-tree.

  boost::dynamic_bitset<> isValidPoint; // Indicates whether a point is valid or has been removed.
  int numValidPoints;                   // The number of non-removed points.
};
