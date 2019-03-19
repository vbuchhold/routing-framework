#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iterator>
#include <limits>
#include <random>
#include <vector>

#include <stocc/stocc.h>
#include <randomc/mersenne.cpp>
#include <randomc/userintf.cpp>
#include <stocc/stoc1.cpp>

#include "DataStructures/Geometry/Point.h"
#include "DataStructures/Geometry/Rectangle.h"
#include "DataStructures/Graph/Graph.h"
#include "Tools/Constants.h"
#include "Tools/OpenMP.h"

// A kd-tree-based approach to choose the closest opportunity with high fitness. In the radiation
// model, each opportunity has a fitness value drawn from some distribution, and each traveler has
// a threshold drawn from the same distribution. Then the traveler chooses the closest opportunity
// with a fitness value higher than their threshold.
template <typename GraphT>
class KDTreeOpportunityChooser {
 public:
  // Constructs an opportunity chooser based on a kd-tree.
  explicit KDTreeOpportunityChooser(const GraphT& graph, const int seed)
      : graph(graph), urbg(seed + omp_get_thread_num() + 1), nrng(seed + omp_get_thread_num()) {
    assert(seed >= 0);
    buildKDTree();
  }

  // Returns the vertex with the closest opportunity with sufficiently high fitness.
  int findClosestOpportunityWithHighFitness(const int src, const int numFitOpportunities) {
    assert(numFitOpportunities > 0);
    source = {graph.latLng(src).longitude(), graph.latLng(src).latitude()};
    boundingBox = {{-INFTY, -INFTY}, {INFTY, INFTY}};
    distToClosestOpportunity = std::numeric_limits<int64_t>::max();
    findClosestOpportunityWithHighFitness(tree.front(), numFitOpportunities);
    assert(boundingBox == Rectangle(Point(-INFTY, -INFTY), Point(INFTY, INFTY)));
    assert(distToClosestOpportunity != std::numeric_limits<int64_t>::max());
    return closestOpportunity;
  }

 private:
  // A (leaf or interior) node in the kd-tree.
  struct Node {
    bool isLeaf;              // A flag that indicates whether this is a leaf or interior node.
    int8_t splitDim;          // The dimension on which we split the record space.
    int32_t splitVal;         // The split value that defines the partition of the record space.
    int32_t leftChild;        // The index of the left child of this node.
    int32_t rightChild;       // The index of the right child of this node.
    int32_t firstRecord;      // The index of the first record in the bucket of this node.
    int32_t lastRecord;       // The index one past the last record in the bucket of this node.
    int32_t numOpportunities; // The number of opportunities in the space represented by this node.
  };

  // A record in the kd-tree, i.e., some vertex with a nonzero number of opportunities.
  struct Record {
    int32_t id;
    Point coordinates;
  };

  static constexpr int BUCKET_SIZE = 16;           // The number of records per bucket.
  static constexpr int RECURSION_THRESHOLD = 1024; // Used to stop the recursion during queries.

  // Builds a kd-tree for all vertices with a nonzero number of opportunities.
  void buildKDTree() {
    FORALL_VERTICES(graph, u)
      if (graph.numOpportunities(u) > 0) {
        Record record;
        record.id = u;
        record.coordinates = {graph.latLng(u).longitude(), graph.latLng(u).latitude()};
        recordsByX.push_back(record);
      }
    Record record;
    record.coordinates = {INFTY, INFTY};
    recordsByX.push_back(record);

    recordsByY.assign(recordsByX.begin(), recordsByX.end());
    tmpStorage.assign(recordsByX.size() - 1, {});
    numOpportunitiesByRecords.assign(recordsByX.size() - 1, 0);
    bucket.assign(recordsByX.size() - 1, 0);

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
      auto numOpportunitiesAtRoot = 0;
      for (auto i = first; i < last; ++i) {
        assert(graph.numOpportunities(recordsByX[i].id) > 0);
        numOpportunitiesByRecords[i] = graph.numOpportunities(recordsByX[i].id);
        numOpportunitiesAtRoot += numOpportunitiesByRecords[i];
      }

      tree[root].isLeaf = true;
      tree[root].splitDim = 0;
      tree[root].splitVal = 0;
      tree[root].leftChild = 0;
      tree[root].rightChild = 0;
      tree[root].firstRecord = first;
      tree[root].lastRecord = last;
      tree[root].numOpportunities = numOpportunitiesAtRoot;
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

    auto numOpportunitiesAtRoot = 0;
    numOpportunitiesAtRoot += tree[tree[root].leftChild].numOpportunities;
    numOpportunitiesAtRoot += tree[tree[root].rightChild].numOpportunities;
    tree[root].numOpportunities = numOpportunitiesAtRoot;
    tree[root].firstRecord = tree[tree[root].leftChild].firstRecord;
    tree[root].lastRecord = tree[tree[root].rightChild].lastRecord;
  }

  // Searches for the closest fit opportunity in the subtree rooted at the specified node.
  void findClosestOpportunityWithHighFitness(const Node& node, const int numFitOpportunities) {
    assert(numFitOpportunities > 0);
    if (node.isLeaf ||
        numFitOpportunities * int64_t{node.lastRecord - node.firstRecord} <= RECURSION_THRESHOLD) {
      handleBaseCase(node, numFitOpportunities);
      return;
    }

    // Draw the number of sufficiently fit opportunities in the left subtree.
    const auto numDraws = numFitOpportunities;
    const auto numSuccesses = tree[node.leftChild].numOpportunities;
    const auto population = node.numOpportunities;
    const auto numFitOpportunitiesOnLeft = nrng.Hypergeometric(numDraws, numSuccesses, population);
    const auto numFitOpportunitiesOnRight = numFitOpportunities - numFitOpportunitiesOnLeft;

    if (source[node.splitDim] <= node.splitVal) {
      // Recursive call on the closer child.
      if (numFitOpportunitiesOnLeft > 0) {
        auto& upperBound = boundingBox.northEast()[node.splitDim];
        const auto tmp = upperBound;
        upperBound = node.splitVal;
        findClosestOpportunityWithHighFitness(tree[node.leftChild], numFitOpportunitiesOnLeft);
        upperBound = tmp;
      }
      // Recursive call on the farther child.
      if (numFitOpportunitiesOnRight > 0) {
        auto& lowerBound = boundingBox.southWest()[node.splitDim];
        const auto tmp = lowerBound;
        lowerBound = node.splitVal;
        if (boundsIntersectBall())
          findClosestOpportunityWithHighFitness(tree[node.rightChild], numFitOpportunitiesOnRight);
        lowerBound = tmp;
      }
    } else {
      // Recursive call on the closer child.
      if (numFitOpportunitiesOnRight > 0) {
        auto& lowerBound = boundingBox.southWest()[node.splitDim];
        const auto tmp = lowerBound;
        lowerBound = node.splitVal;
        findClosestOpportunityWithHighFitness(tree[node.rightChild], numFitOpportunitiesOnRight);
        lowerBound = tmp;
      }
      // Recursive call on the farther child.
      if (numFitOpportunitiesOnLeft > 0) {
        auto& upperBound = boundingBox.northEast()[node.splitDim];
        const auto tmp = upperBound;
        upperBound = node.splitVal;
        if (boundsIntersectBall())
          findClosestOpportunityWithHighFitness(tree[node.leftChild], numFitOpportunitiesOnLeft);
        upperBound = tmp;
      }
    }
  }

  // Samples numFitOpportunities opportunities from the record space represented by the given node
  // and checks for each if it improves the closest fit opportunity so far encountered.
  void handleBaseCase(const Node& node, const int numFitOpportunities) {
    assert(numFitOpportunities <= node.numOpportunities);
    const auto beginning = numOpportunitiesByRecords.begin();
    std::copy(beginning + node.firstRecord, beginning + node.lastRecord, bucket.begin());
    for (auto i = 1; i <= numFitOpportunities; ++i) {
      const auto rank = std::uniform_int_distribution<>(0, node.numOpportunities - i)(urbg);
      auto j = 0;
      auto partialSum = bucket.front();
      while (partialSum <= rank) partialSum += bucket[++j];
      assert(j < node.lastRecord - node.firstRecord);
      assert(partialSum <= node.numOpportunities);
      assert(bucket[j] > 0);
      --bucket[j];
      const auto& record = buckets[node.firstRecord + j];
      const auto dist = source.getSquaredEuclideanDistanceTo(record.coordinates);
      if (dist < distToClosestOpportunity) {
        closestOpportunity = record.id;
        distToClosestOpportunity = dist;
      }
    }
  }

  // Returns true if the current bounding box intersects the ball centered at the source whose
  // radius is equal to the distance to the closest fit opportunity so far encountered.
  bool boundsIntersectBall() const noexcept {
    int64_t distToBox = 0;
    if (source.x() < boundingBox.southWest().x()) {
      const int64_t diff = source.x() - boundingBox.southWest().x();
      distToBox += diff * diff;
    } else if (source.x() > boundingBox.northEast().x()) {
      const int64_t diff = source.x() - boundingBox.northEast().x();
      distToBox += diff * diff;
    }
    if (source.y() < boundingBox.southWest().y()) {
      const int64_t diff = source.y() - boundingBox.southWest().y();
      distToBox += diff * diff;
    } else if (source.y() > boundingBox.northEast().y()) {
      const int64_t diff = source.y() - boundingBox.northEast().y();
      distToBox += diff * diff;
    }
    return distToBox < distToClosestOpportunity;
  }

  const GraphT& graph;   // The network we work on.
  std::minstd_rand urbg; // A uniform random bit generator.
  StochasticLib1 nrng;   // A nonuniform random number generator.

  std::vector<Node> tree;                         // The nodes in the kd-tree.
  std::vector<Record> buckets;                    // The buckets of the leaves concatenated.
  std::vector<int32_t> numOpportunitiesByRecords; // The number of opportunities for each record.

  Point source;                     // The source for which we choose the corresponding target.
  Rectangle boundingBox;            // The bounds of the subspace represented by the current node.
  std::vector<int32_t> bucket;      // The number of opportunities for records in current bucket.
  int closestOpportunity;           // The closest fit opportunity so far encountered.
  int64_t distToClosestOpportunity; // The distance to the closest fit opportunity.

  std::vector<Record> recordsByX; // During construction: records ordered by the x-coordinate.
  std::vector<Record> recordsByY; // During construction: records ordered by the y-coordinate.
  std::vector<Record> tmpStorage; // During construction: records larger than the split value.
};
