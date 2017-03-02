#pragma once

namespace kdtree {

// A result set for a circular range query in a kd-tree. For each point within the query range, it
// invokes a user-provided hook function.
template <typename PointSetAdapterT, typename DistanceT, typename PointFoundT>
class CircularRangeResultSet {
 public:
  // Constructs a result set for a circular range query with the specified radius.
  CircularRangeResultSet(const PointSetAdapterT& pointSet, const DistanceT radius,
                         PointFoundT pointFound)
      : pointSet(pointSet), radius(radius), pointFound(pointFound) {}

  // Returns the maximum distance a point must have to be potentially inserted into the result.
  DistanceT worstDist() const {
    return radius + 1; // +1 since a point's distance should be less than OR EQUAL TO the radius.
  }

  // In the case of k-nearest-neighbor queries, returns true if k points were found.
  bool full() const {
    return true;
  }

  // Invoked for each point within a distance of at most worstDist() from the query point.
  void addPoint(const DistanceT dist, const int idx) const {
    pointFound(pointSet.getIndexInPointSet(idx), dist);
  }

 private:
  const PointSetAdapterT& pointSet; // An adapter to the set of points that should be searched.
  const DistanceT radius;           // The radius of the query range.
  PointFoundT pointFound;           // A hook invoked for each point within the query range.
};

}
