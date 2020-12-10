#pragma once

#include <algorithm>
#include <cstdint>
#include <utility>
#include <vector>

#include "Algorithms/Buckets/BucketEntry.h"
#include "Algorithms/Buckets/StaticBucketContainer.h"

namespace buckets {

// This class represents an entry in the bucket of a vertex v that stores the hub v explicitly.
struct BucketEntryWithHub : public BucketEntry {
  // Creates an entry that represents a path between hub and targetId of length distToTarget.
  BucketEntryWithHub(const int hub, const int targetId, const int distToTarget) noexcept
      : BucketEntry(targetId, distToTarget), hub(hub) {}

  int32_t hub; // The hub whose bucket contains this entry.
};

// An (abused) pruning criterion for the reverse CH search that always returns false but creates
// and collects bucket entries for all vertices in the search space as a side effect.
struct CollectSearchSpace {
  // Creates an instance of this pruning criterion.
  CollectSearchSpace(const int& targetId, std::vector<BucketEntryWithHub>& searchSpace) noexcept
      : targetId(targetId), searchSpace(searchSpace) {}

  // Always returns false but creates and collects bucket entries as a side effect.
  template <typename DistLabelT, typename DistLabelContainerT>
  bool operator()(const int v, DistLabelT& distToV, const DistLabelContainerT& /*distLabels*/) {
    searchSpace.emplace_back(v, targetId, distToV[0]);
    return false;
  }

  const int& targetId;                          // The target for which we are generating entries.
  std::vector<BucketEntryWithHub>& searchSpace; // A dynamic array that collects the entries.
};

}

// This class implements the construction of static buckets for a given set of targets. Note that
// the class cannot be instantiated because the details of the CH search used depend on whether we
// use standard or customizable CH. See the derived classes BCHBucketBuilder and BCCHBucketBuilder
// for concrete bucket builders.
template <typename CHSearchT>
class BucketBuilder {
 public:
  // Creates a bucket builder that uses the specified CH search to create and collect the entries.
  BucketBuilder(const int numVertices, CHSearchT&& search)
      : numVertices(numVertices), search(std::move(search)) {}

  // Returns buckets that store the reverse CH search spaces of the specified set of targets.
  template <typename IteratorT>
  StaticBucketContainer buildBucketsFor(
      const IteratorT firstTarget, const IteratorT lastTarget, const bool sorted = false) {
    targetId = 0;
    searchSpaces.clear();
    for (auto t = firstTarget; t != lastTarget; ++t) {
      search.run(*t);
      ++targetId;
    }

    if (sorted)
      std::sort(searchSpaces.begin(), searchSpaces.end(), [&](const auto& lhs, const auto& rhs) {
        return lhs.hub < rhs.hub || (lhs.hub == rhs.hub && lhs.distToTarget < rhs.distToTarget);
      });
    else
      std::sort(searchSpaces.begin(), searchSpaces.end(), [&](const auto& lhs, const auto& rhs) {
        return lhs.hub < rhs.hub;
      });

    StaticBucketContainer cont;
    cont.bucketPositions.assign(numVertices + 1, 0);
    cont.entries.assign(searchSpaces.begin(), searchSpaces.end());
    searchSpaces.emplace_back(numVertices, 0, 0);
    for (auto v = 0, i = 0; v < cont.bucketPositions.size(); ++v) {
      while (searchSpaces[i].hub < v) ++i;
      cont.bucketPositions[v] = i;
    }
    return cont;
  }

 protected:
  using CHSearch = CHSearchT;
  using SearchSpaceContainer = std::vector<buckets::BucketEntryWithHub>;

  int targetId;                      // The target for which we are generating bucket entries.
  SearchSpaceContainer searchSpaces; // A dynamic array that collects the bucket entries.

 private:
  const int numVertices; // The number of vertices in the CH.
  CHSearch search;       // The CH search that creates and collects the bucket entries.
};
