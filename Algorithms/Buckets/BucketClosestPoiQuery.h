#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

#include "Algorithms/Buckets/StaticBucketContainer.h"
#include "DataStructures/Queues/AddressableKHeap.h"
#include "Tools/Constants.h"

namespace buckets {

// The stopping criterion for the forward search. We can stop when the search scans a vertex whose
// distance label is no smaller than the distance to the k-th closest POI seen so far.
struct StopWhenFrontierFartherThanClosestPois {
  // Creates an instance of the stopping criterion for the forward CH search.
  StopWhenFrontierFartherThanClosestPois(const AddressableQuadheap& currentlyClosestPois) noexcept
      : currentlyClosestPois(currentlyClosestPois) {}

  // Returns true if the search can be stopped at v.
  template <typename DistLabelT, typename DistLabelContT>
  bool operator()(const int /*v*/, DistLabelT& distToV, const DistLabelContT& /*distLabels*/) {
    return distToV[0] >= -currentlyClosestPois.minKey();
  }

  const AddressableQuadheap& currentlyClosestPois; // The closest POIs encountered so far.
};

// An (abused) pruning criterion that always returns false but scans the buckets as a side effect.
struct ScanBucket {
  // Creates an instance of this pruning criterion.
  ScanBucket(
      const int& numPoisToBeReported, const StaticBucketContainer* const& cont,
      AddressableQuadheap& currentlyClosestPois) noexcept
      : numPoisToBeReported(numPoisToBeReported),
        bucketContainer(cont),
        currentlyClosestPois(currentlyClosestPois) {}

  // Always returns false but scans the bucket of v as a side effect.
  template <typename DistLabelT, typename DistLabelContainerT>
  bool operator()(const int v, DistLabelT& distToV, const DistLabelContainerT& /*distLabels*/) {
    assert(bucketContainer != nullptr);
    for (const auto& entry : bucketContainer->getBucketOf(v)) {
      const auto distViaV = distToV[0] + entry.distToTarget;
      if (distViaV >= -currentlyClosestPois.minKey())
        break;
      if (currentlyClosestPois.contains(entry.targetId)) {
        currentlyClosestPois.increaseKeyIfPossible(entry.targetId, -distViaV);
      } else {
        if (currentlyClosestPois.size() == numPoisToBeReported) {
          int poi, distToPoi;
          currentlyClosestPois.deleteMin(poi, distToPoi);
        }
        currentlyClosestPois.insert(entry.targetId, -distViaV);
      }
      assert(currentlyClosestPois.size() <= numPoisToBeReported);
    }
    return false;
  }

  const int& numPoisToBeReported;                      // The number of POIs to be reported.
  const StaticBucketContainer* const& bucketContainer; // The reverse search spaces of the POIs.
  AddressableQuadheap& currentlyClosestPois;           // The closest POIs encountered so far.
};

}

// This class implements a bucket-based closest-POI query in a road network. It works in two phases.
// Given a set of POIs, the caller must first build a POI index (static buckets for the set of POIs)
// by invoking the method buildPoiIndexFor(). The index is then used to run closest-POI queries by
// invoking the method findClosestPois().
//
// Note that the class cannot be instantiated because the details of the CH searches that build and
// scan the buckets depend on whether we use standard or customizable CH. See the derived classes
// BCHClosestPoiQuery and BCCHClosestPoiQuery for concrete queries.
template <typename BucketBuilderT, typename ForwardCHSearchT>
class BucketClosestPoiQuery {
 public:
  // Precomputed auxiliary information to accelerate POI queries.
  struct PoiIndex {
    friend class BucketClosestPoiQuery;

   private:
    // Builds the POI index for the specified set of POI vertices.
    PoiIndex(const std::vector<int32_t>& pointsOfInterest, StaticBucketContainer&& bucketCont)
        : pointsOfInterest(pointsOfInterest), bucketContainer(std::move(bucketCont)) {}

    const std::vector<int32_t>& pointsOfInterest; // The POI vertices in increasing order of ID.
    const StaticBucketContainer bucketContainer;  // The reverse search spaces of the POI vertices.
  };

  // A point of interest returned by a POI query.
  struct Poi {
    // Returns true if this POI is closer to the source than the specified POI.
    constexpr bool operator<(const Poi& rhs) noexcept {
      return dist < rhs.dist;
    }

    int vertex; // The vertex that contains this POI.
    int dist;   // The shortest-path distance from the source to this POI.
  };

  // Creates a bucket-based closest-POI query that uses the given bucket builder and forward search.
  BucketClosestPoiQuery(
      const int numVertices, BucketBuilderT* const builder, ForwardCHSearchT&& search)
      : bucketContainer(nullptr),
        currentlyClosestPois(numVertices + 1),
        bucketBuilder(builder),
        search(std::move(search)) {}

  // Builds the POI index for the specified set of POI vertices.
  PoiIndex buildPoiIndexFor(const std::vector<int32_t>& pointsOfInterest) {
    return {pointsOfInterest,
            bucketBuilder->buildBucketsFor(pointsOfInterest.begin(), pointsOfInterest.end(), true)};
  }

  // Returns the k closest POI vertices to s.
  const std::vector<Poi>& findClosestPois(const int s, const PoiIndex& idx, const int k = 1) {
    numPoisToBeReported = k;
    bucketContainer = &idx.bucketContainer;
    currentlyClosestPois.insert(idx.pointsOfInterest.size(), -INFTY);
    search.run(s);

    if (currentlyClosestPois.minId() == idx.pointsOfInterest.size()) {
      int poi, distToPoi;
      currentlyClosestPois.deleteMin(poi, distToPoi);
    }
    if (currentlyClosestPois.size() < closestPois.size())
      closestPois.erase(closestPois.begin() + currentlyClosestPois.size(), closestPois.end());
    else
      closestPois.insert(closestPois.end(), currentlyClosestPois.size() - closestPois.size(), {});
    for (auto iter = closestPois.rbegin(); iter < closestPois.rend(); ++iter) {
      int poi, distToPoi;
      currentlyClosestPois.deleteMin(poi, distToPoi);
      iter->vertex = idx.pointsOfInterest[poi];
      iter->dist = -distToPoi;
    }
    assert(currentlyClosestPois.empty());
    assert(std::is_sorted(closestPois.begin(), closestPois.end()));
    return closestPois;
  }

 protected:
  using BucketBuilder = BucketBuilderT;
  using ForwardCHSearch = ForwardCHSearchT;

  int numPoisToBeReported;                      // The number of POIs to be reported.
  const StaticBucketContainer* bucketContainer; // The reverse search spaces of the set of POIs.
  AddressableQuadheap currentlyClosestPois;     // The closest POIs, with the farthest one on top.

 private:
  const std::unique_ptr<BucketBuilder> bucketBuilder; // The facility for building the buckets.
  ForwardCHSearch search;                             // The CH search that scans the buckets.
  std::vector<Poi> closestPois;                       // The result of the last query, close to far.
};
