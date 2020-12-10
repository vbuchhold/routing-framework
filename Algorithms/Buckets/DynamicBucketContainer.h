#pragma once

#include <algorithm>
#include <cassert>
#include <vector>

#include "DataStructures/Utilities/DynamicRagged2DArrays.h"

// This class maintains a dynamic bucket for each vertex for bucket-based CH searches. We store all
// bucket entries in a single dynamic value array. The entries in the same bucket are stored
// consecutively in memory, in no particular order. In addition, an index array stores the starting
// and ending point of each bucket's value block in the value array.
//
// When we remove an entry from a bucket, we fill the resulting hole in the value array with the
// rightmost value in the bucket's value block, and decrement the block's ending point in the index
// array. Consider an insertion of an entry into a bucket. If the element immediately after the
// bucket's value block is a hole, the new entry fills this hole. Analogously, if the element before
// the value block is a hole, the new entry fills that hole. Otherwise, we move the entire value
// block to the end of the value array, and additionally insert a number of holes after the value
// block (the number is a constant fraction of the block size). Then, there is a hole after the
// block, and we proceed as described above.
template <typename BucketEntryT>
class DynamicBucketContainer {
 public:
  using Bucket = ConstantValueBlock<BucketEntryT>;

  // Constructs a container that can maintain buckets for the specified number of vertices.
  explicit DynamicBucketContainer(const int numVertices) : bucketPositions(numVertices) {
    assert(numVertices >= 0);
  }

  // Returns the bucket of the specified vertex.
  Bucket getBucketOf(const int vertex) const {
    assert(vertex >= 0); assert(vertex < bucketPositions.size());
    const auto& pos = bucketPositions[vertex];
    return Bucket(entries.begin() + pos.start, entries.begin() + pos.end);
  }

  // Inserts the given entry into the bucket of the specified vertex.
  bool insert(const int vertex, const BucketEntryT& entry) {
    insertion(vertex, entry, bucketPositions, entries);
    return true;
  }

  // Removes the entry for targetId from the bucket of the specified vertex.
  bool remove(const int vertex, const int targetId) {
    assert(vertex >= 0); assert(vertex < bucketPositions.size());
    numEntriesVisited = 0;
    const auto& pos = bucketPositions[vertex];
    for (auto i = pos.start; i < pos.end; ++i) {
      ++numEntriesVisited;
      if (entries[i].targetId == targetId) {
        removal(vertex, i - pos.start, bucketPositions, entries);
        return true;
      }
    }
    return false;
  }

  // Removes all entries from all buckets.
  void clear() {
    for (auto& bucketPos : bucketPositions)
      bucketPos.end = bucketPos.start;
    std::fill(entries.begin(), entries.end(), BucketEntryT());
  }

  // Returns the number of bucket entries visited during the last remove operation.
  int getNumEntriesVisited() const noexcept {
    return numEntriesVisited;
  }

 private:
  using BucketPosition = ValueBlockPosition;

  std::vector<BucketPosition> bucketPositions;
  std::vector<BucketEntryT> entries;
  int numEntriesVisited;
};
