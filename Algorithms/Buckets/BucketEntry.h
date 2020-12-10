#pragma once

#include <cstdint>
#include <limits>

// This class represents an entry in the bucket of a vertex v. It can be thought of as a shortcut
// from v to targetId with length distToTarget. Note that the terminology was developed with
// one-to-many queries in mind. In many-to-one queries, targetId stores a source ID and distToTarget
// stores the distance from the corresponding source to v.
struct BucketEntry {
  BucketEntry() noexcept = default;

  BucketEntry(const int targetId, const int distToTarget) noexcept
      : targetId(targetId), distToTarget(distToTarget) {}

  constexpr bool operator==(const BucketEntry& rhs) const noexcept {
    return targetId == rhs.targetId;
  }

  int32_t targetId = std::numeric_limits<int32_t>::max();
  int32_t distToTarget;
};
