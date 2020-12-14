#pragma once

#include <cassert>
#include <cstdint>
#include <vector>

#include "Algorithms/Buckets/BucketEntry.h"
#include "DataStructures/Utilities/DynamicRagged2DArrays.h"

// This class maintains a static bucket for each vertex for bucket-based CH searches. We store all
// bucket entries in a single value array. The entries in the same bucket are stored consecutively
// in memory. In addition, an index array stores the starting position of each bucket's value block
// in the value array.
class StaticBucketContainer {
  template <typename>
  friend class BucketBuilder;

 public:
  using Bucket = ConstantValueBlock<BucketEntry>;

  // Returns the bucket of v.
  Bucket getBucketOf(const int v) const {
    assert(v >= 0); assert(v < bucketPositions.size() - 1);
    return {entries.begin() + bucketPositions[v], entries.begin() + bucketPositions[v + 1]};
  }

  // Returns the space (in bytes) consumed by these buckets.
  int spaceConsumption() const noexcept {
    return bucketPositions.size() * sizeof(int32_t) + entries.size() * sizeof(BucketEntry);
  }

 private:
  // Default constructor.
  StaticBucketContainer() noexcept = default;

  std::vector<int32_t> bucketPositions; // The starting position of each bucket in the entry array.
  std::vector<BucketEntry> entries;     // With entries in the same bucket stored consecutively.
};
