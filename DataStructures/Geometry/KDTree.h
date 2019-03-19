#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iterator>
#include <limits>
#include <vector>

#include "DataStructures/Geometry/Point.h"
#include "DataStructures/Geometry/Rectangle.h"
#include "Tools/Constants.h"

// A kd-tree for the problem of finding nearest neighbors.
class KDTree {
 public:
  // Constructs a kd-tree that stores the specified set of points.
  explicit KDTree(const std::vector<Point>& points) {
    buildKDTree(points);
  }

  // Returns the point that is closest to the query point while being no farther than maxDist.
  int findClosestPoint(const Point& query, int32_t maxDist = std::numeric_limits<int32_t>::max()) {
    assert(maxDist >= 0);
    queryPoint = query;
    boundingBox = {{-INFTY, -INFTY}, {INFTY, INFTY}};
    closestPoint = INVALID_ID;
    distToClosestPoint = int64_t{maxDist} * maxDist + 1;
    findClosestPoint(tree.front());
    assert(boundingBox == Rectangle(Point(-INFTY, -INFTY), Point(INFTY, INFTY)));
    return closestPoint;
  }

 private:
  // A (leaf or interior) node in the kd-tree.
  struct Node {
    bool isLeaf;           // A flag that indicates whether this is a leaf or interior node.
    int8_t splitDim;       // The dimension on which we split the record space.
    int32_t splitVal;      // The split value that defines the partition of the record space.
    union {
      int32_t leftChild;   // The index of the left child of this node.
      int32_t firstRecord; // The index of the first record in the bucket of this node.
    };
    union {
      int32_t rightChild;  // The index of the right child of this node.
      int32_t lastRecord;  // The index one past the last record in the bucket of this node.
    };
  };

  // A record in the kd-tree, i.e., a point together with its ID.
  struct Record {
    int32_t id;
    Point coordinates;
  };

  static constexpr int BUCKET_SIZE = 16; // The number of records per bucket.

  // Builds a kd-tree that stores the specified set of points.
  void buildKDTree(const std::vector<Point>& points) {
    recordsByX.assign(points.size() + 1, {});
    for (auto i = 0; i < points.size(); ++i) {
      recordsByX[i].id = i;
      recordsByX[i].coordinates = points[i];
    }
    recordsByX.back().coordinates = {INFTY, INFTY};
    recordsByY.assign(recordsByX.begin(), recordsByX.end());
    tmpStorage.assign(points.size(), {});

    std::stable_sort(recordsByX.begin(), recordsByX.end() - 1, [](const auto& a, const auto& b) {
      return a.coordinates.x() < b.coordinates.x();
    });
    std::stable_sort(recordsByY.begin(), recordsByY.end() - 1, [](const auto& a, const auto& b) {
      return a.coordinates.y() < b.coordinates.y();
    });

    buildKDTree(0, recordsByX.size() - 1, 0);
    recordsByX.swap(buckets);
    recordsByY = std::vector<Record>();
    tmpStorage = std::vector<Record>();
    buckets.pop_back();
  }

  // Builds a kd-tree for the records in the range [first, last).
  void buildKDTree(const int first, const int last, const int splitDim) {
    assert(0 <= first); assert(first < last); assert(last < recordsByX.size());
    assert(splitDim == 0 || splitDim == 1);
    const auto root = tree.size();
    tree.emplace_back();

    if (last - first <= BUCKET_SIZE) {
      tree[root].isLeaf = true;
      tree[root].splitDim = 0;
      tree[root].splitVal = 0;
      tree[root].firstRecord = first;
      tree[root].lastRecord = last;
      return;
    }

    tree[root].isLeaf = false;
    tree[root].splitDim = splitDim;

    if (splitDim == 0) {
      // Determine the split value.
      auto middle = (first + last) / 2;
      tree[root].splitVal = recordsByX[middle].coordinates[splitDim];
      while (recordsByX[middle].coordinates[splitDim] == tree[root].splitVal) ++middle;
      assert(middle != last);

      // Partition the record list ordered by the y-coordinate with respect to the x-coordinate.
      auto smallRecords = recordsByY.begin() + first;
      auto largeRecords = tmpStorage.begin();
      for (auto iter = smallRecords; iter != recordsByY.begin() + last; ++iter)
        if (iter->coordinates[splitDim] <= tree[root].splitVal)
          *smallRecords++ = *iter;
        else
          *largeRecords++ = *iter;
      assert(std::distance(recordsByY.begin(), smallRecords) == middle);
      assert(std::distance(tmpStorage.begin(), largeRecords) == last - middle);
      std::copy(tmpStorage.begin(), largeRecords, smallRecords);

      // Recurse on the two subproblems.
      tree[root].leftChild = tree.size();
      buildKDTree(first, middle, !splitDim);
      tree[root].rightChild = tree.size();
      buildKDTree(middle, last, !splitDim);
    } else {
      // Determine the split value.
      auto middle = (first + last) / 2;
      tree[root].splitVal = recordsByY[middle].coordinates[splitDim];
      while (recordsByY[middle].coordinates[splitDim] == tree[root].splitVal) ++middle;
      assert(middle != last);

      // Partition the record list ordered by the x-coordinate with respect to the y-coordinate.
      auto smallRecords = recordsByX.begin() + first;
      auto largeRecords = tmpStorage.begin();
      for (auto iter = smallRecords; iter != recordsByX.begin() + last; ++iter)
        if (iter->coordinates[splitDim] <= tree[root].splitVal)
          *smallRecords++ = *iter;
        else
          *largeRecords++ = *iter;
      assert(std::distance(recordsByX.begin(), smallRecords) == middle);
      assert(std::distance(tmpStorage.begin(), largeRecords) == last - middle);
      std::copy(tmpStorage.begin(), largeRecords, smallRecords);

      // Recurse on the two subproblems.
      tree[root].leftChild = tree.size();
      buildKDTree(first, middle, !splitDim);
      tree[root].rightChild = tree.size();
      buildKDTree(middle, last, !splitDim);
    }
  }

  // Searches for the closest point in the subtree rooted at the specified node.
  void findClosestPoint(const Node& node) {
    if (node.isLeaf) {
      handleBaseCase(node);
      return;
    }

    if (queryPoint[node.splitDim] <= node.splitVal) {
      // Recursive call on the closer child.
      auto& upperBound = boundingBox.northEast()[node.splitDim];
      auto tmp = upperBound;
      upperBound = node.splitVal;
      findClosestPoint(tree[node.leftChild]);
      upperBound = tmp;

      // Recursive call on the farther child.
      auto& lowerBound = boundingBox.southWest()[node.splitDim];
      tmp = lowerBound;
      lowerBound = node.splitVal;
      if (boundsIntersectBall())
        findClosestPoint(tree[node.rightChild]);
      lowerBound = tmp;
    } else {
      // Recursive call on the closer child.
      auto& lowerBound = boundingBox.southWest()[node.splitDim];
      auto tmp = lowerBound;
      lowerBound = node.splitVal;
      findClosestPoint(tree[node.rightChild]);
      lowerBound = tmp;

      // Recursive call on the farther child.
      auto& upperBound = boundingBox.northEast()[node.splitDim];
      tmp = upperBound;
      upperBound = node.splitVal;
      if (boundsIntersectBall())
        findClosestPoint(tree[node.leftChild]);
      upperBound = tmp;
    }
  }

  // Checks for each point in the record space represented by the specified node if it improves the
  // closest point so far encountered.
  void handleBaseCase(const Node& node) {
    for (auto i = node.firstRecord; i < node.lastRecord; ++i) {
      const auto dist = queryPoint.getSquaredEuclideanDistanceTo(buckets[i].coordinates);
      if (dist < distToClosestPoint) {
        closestPoint = buckets[i].id;
        distToClosestPoint = dist;
      }
    }
  }

  // Returns true if the current bounding box intersects the ball centered at the query point whose
  // radius is equal to the distance to the closest point so far encountered.
  bool boundsIntersectBall() const noexcept {
    int64_t distToBox = 0;
    if (queryPoint.x() < boundingBox.southWest().x()) {
      const int64_t diff = queryPoint.x() - boundingBox.southWest().x();
      distToBox += diff * diff;
    } else if (queryPoint.x() > boundingBox.northEast().x()) {
      const int64_t diff = queryPoint.x() - boundingBox.northEast().x();
      distToBox += diff * diff;
    }
    if (queryPoint.y() < boundingBox.southWest().y()) {
      const int64_t diff = queryPoint.y() - boundingBox.southWest().y();
      distToBox += diff * diff;
    } else if (queryPoint.y() > boundingBox.northEast().y()) {
      const int64_t diff = queryPoint.y() - boundingBox.northEast().y();
      distToBox += diff * diff;
    }
    return distToBox < distToClosestPoint;
  }

  std::vector<Node> tree;      // The nodes in the kd-tree.
  std::vector<Record> buckets; // The buckets of the leaves concatenated.

  Point queryPoint;           // The query point for which we search the nearest neighbor.
  Rectangle boundingBox;      // The bounds of the record space represented by the current node.
  int closestPoint;           // The closest point so far encountered.
  int64_t distToClosestPoint; // The distance to the closest point.

  std::vector<Record> recordsByX; // During construction: records ordered by the x-coordinate.
  std::vector<Record> recordsByY; // During construction: records ordered by the y-coordinate.
  std::vector<Record> tmpStorage; // During construction: records larger than the split value.
};
