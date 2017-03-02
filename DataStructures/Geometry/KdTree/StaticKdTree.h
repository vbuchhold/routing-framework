#pragma once

#include <cassert>

#include <nanoflann.hpp>

#include "DataStructures/Geometry/KdTree/Adapters/PointSetAdapter.h"
#include "DataStructures/Geometry/KdTree/Metrics/EuclideanMetric.h"
#include "DataStructures/Geometry/KdTree/ResultSets/CircularRangeResultSet.h"
#include "DataStructures/Geometry/Point.h"

// A static kd-tree for sets of two-dimensional points. It supports efficient circular range queries
// and k-nearest-neighbor finding. The metric used to compute distances (Euclidean or Manhattan) is
// specified by a template parameter.
template <
    typename PointSetAdapterT = kdtree::PointSetAdapter,
    template <typename> class MetricT = kdtree::EuclideanMetric>
class StaticKdTree {
 private:
  using Metric = MetricT<PointSetAdapterT>;       // The metric used to compute distances.
  using Distance = typename Metric::DistanceType; // The distance type.

 public:
  // Constructs a static kd-tree using the specified PointSet adapter.
  StaticKdTree(const PointSetAdapterT& points) : pointSet(points), kdTree(2, pointSet) {
    rebuild();
  }

  // Rebuilds this kd-tree.
  void rebuild() {
    kdTree.buildIndex();
  }

  // Invokes for each point within the specified circular query range the given hook function.
  template <typename PointFoundT>
  void circularRangeQuery(const Point& p, const Distance radius, PointFoundT pointFound) const {
    using ResultSet = kdtree::CircularRangeResultSet<PointSetAdapterT, Distance, PointFoundT>;
    ResultSet res(pointSet, radius, pointFound);
    const int queryPoint[] = {p.getX(), p.getY()};
    kdTree.findNeighbors(res, queryPoint, nanoflann::SearchParams());
  }

  // Locates the k points closest to the query point p.
  void nearestNeighborQuery(const Point& p, const int k, int* indices, Distance* distances) const {
    assert(k > 0);
    const int queryPoint[] = {p.getX(), p.getY()};
    const int numFound = kdTree.knnSearch(queryPoint, k, indices, distances);
    for (int i = 0; i < numFound; ++i)
      indices[i] = pointSet.getIndexInPointSet(indices[i]);
  }

 private:
  using KdTree = nanoflann::KDTreeSingleIndexAdaptor<Metric, PointSetAdapterT, 2, int>;

  const PointSetAdapterT& pointSet; // An adapter to the set of points that should be searched.
  KdTree kdTree;                    // The underlying nanoflann kd-tree.
};
