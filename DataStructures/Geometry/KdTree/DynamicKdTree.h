#pragma once

#include <cassert>

#include <boost/dynamic_bitset.hpp>

#include "DataStructures/Geometry/KdTree/Metrics/EuclideanMetric.h"
#include "DataStructures/Geometry/KdTree/StaticKdTree.h"

// A dynamic kd-tree for sets of two-dimensional points. It supports efficient circular range
// queries and k-nearest-neighbor finding. The metric used to compute distances (Euclidean or
// Manhattan) is specified by a template parameter.
template <
    typename PointSetAdapterT = kdtree::PointSetAdapter,
    template <typename> class MetricT = kdtree::EuclideanMetric>
class DynamicKdTree {
 private:
  using Distance = typename MetricT<PointSetAdapterT>::DistanceType; // The distance type.

 public:
  // Constructs a dynamic kd-tree using the specified PointSet adapter.
  DynamicKdTree(const PointSetAdapterT& initialPointSet)
      : initialPointSet(initialPointSet),
        currentPointSet(currentPoints),
        currentKdTree(currentPointSet) {
    rebuild();
  }

  // Rebuilds this kd-tree.
  void rebuild() {
    numValidPoints = initialPointSet.kdtree_get_point_count();
    isValidPoint.clear();
    isValidPoint.resize(numValidPoints, true);

    currentPoints.clear();
    currentPoints.reserve(numValidPoints);
    indexInInitialSet.resize(numValidPoints);
    for (int idx = 0; idx < numValidPoints; ++idx) {
      currentPoints.insert(initialPointSet[idx]);
      indexInInitialSet[idx] = idx;
    }
    currentKdTree.rebuild();
  }

  // Returns true if this kd-tree contains the point with index idx.
  bool contains(const int idx) const {
    assert(idx >= 0); assert(idx < isValidPoint.size());
    return isValidPoint[idx];
  }

  // Invokes for each point within the specified circular query range the given hook function.
  template <typename PointFoundT>
  void circularRangeQuery(const Point& p, const Distance radius, PointFoundT pointFound) const {
    currentKdTree.circularRangeQuery(p, radius,
                                     [this, pointFound](const int idx, const Distance dist) {
      assert(idx >= 0); assert(idx < indexInInitialSet.size());
      const int idxInInitialSet = indexInInitialSet[idx];
      if (isValidPoint[idxInInitialSet])
        pointFound(idxInInitialSet, dist);
    });
  }

  // Removes the point with index idx from this kd-tree.
  void remove(const int idx) {
    assert(idx >= 0); assert(idx < isValidPoint.size());
    assert(isValidPoint[idx]);
    isValidPoint[idx] = false;
    --numValidPoints;

    // Whenever the number of valid points becomes half the size of the kd-tree, we rebuild it.
    if (2 * numValidPoints <= currentPointSet.kdtree_get_point_count()) {
      currentPoints.clear();
      indexInInitialSet.resize(numValidPoints);
      for (int idx = isValidPoint.find_first(); idx < boost::dynamic_bitset<>::npos;
           idx = isValidPoint.find_next(idx)) {
        assert(currentPoints.size() < numValidPoints);
        currentPoints.insert(initialPointSet[idx]);
        indexInInitialSet[currentPoints.size() - 1] = idx;
      }
      assert(currentPoints.size() == numValidPoints);
      currentKdTree.rebuild();
    }
  }

 private:
  using StaticKdTreeT = StaticKdTree<kdtree::PointSetAdapter, MetricT>;

  const PointSetAdapterT& initialPointSet;       // An adapter to the initial set of points.
  PointSet currentPoints;                        // The points currently stored in the kd-tree.
  const kdtree::PointSetAdapter currentPointSet; // An adapter to the current set of points.
  StaticKdTreeT currentKdTree;                   // A static kd-tree for the current set of points.
  std::vector<int> indexInInitialSet;            // For each point p, p's index in the initial set.

  boost::dynamic_bitset<> isValidPoint; // Indicates whether a point is valid or has been removed.
  int numValidPoints;                   // The number of non-removed points.
};
