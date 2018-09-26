#pragma once

#include <cassert>
#include <cstdint>
#include <numeric>
#include <vector>

#include <omp.h>

#include "DataStructures/Containers/BitVector.h"
#include "Tools/Bitwise.h"

// This class maps a set of n IDs in the range 0..m - 1 to the range 0..n - 1, preserving the
// relative order. We call the IDs in the larger range global IDs, and the IDs in the smaller
// range local IDs. Note that the local IDs are sequential.
template <int k = 64>
class LocalIdMap {
 public:
  // Constructs a map from the specified set of n global IDs to the range 0..n - 1.
  explicit LocalIdMap(const BitVector& mappedGlobalIds)
      : mappedGlobalIds(mappedGlobalIds),
        numMappedGlobalIds((mappedGlobalIds.size() + k - 1) / k + 1) {
    constexpr int blocksPerRank = k / BitVector::BITS_PER_BLOCK;
    std::vector<int> sums;

    #pragma omp parallel
    {
      #pragma omp single
      sums.assign(omp_get_num_threads() + 1, 0);

      int sum = 0;

      #pragma omp for schedule(static) nowait
      for (int i = 0; i < (mappedGlobalIds.size() + BitVector::BITS_PER_BLOCK - 1) / k; ++i) {
        for (int j = 0; j < blocksPerRank; ++j)
          sum += bitCount(mappedGlobalIds.block(i * blocksPerRank + j));
        numMappedGlobalIds[i + 1] = sum;
      }

      if (omp_get_num_threads() > 1) {
        sums[omp_get_thread_num() + 1] = sum;

        #pragma omp barrier

        #pragma omp single
        std::partial_sum(sums.begin(), sums.end(), sums.begin());

        #pragma omp for schedule(static) nowait
        for (int i = 0; i < (mappedGlobalIds.size() + BitVector::BITS_PER_BLOCK - 1) / k; ++i)
          numMappedGlobalIds[i + 1] += sums[omp_get_thread_num()];
      }
    }

    const int i = (mappedGlobalIds.size() + BitVector::BITS_PER_BLOCK - 1) / k;
    if (i + 1 < numMappedGlobalIds.size()) {
      numMappedGlobalIds[i + 1] = numMappedGlobalIds[i];
      for (int j = i * blocksPerRank; j < mappedGlobalIds.numBlocks(); ++j)
        numMappedGlobalIds[i + 1] += bitCount(mappedGlobalIds.block(j));
    }
  }

  // Returns the size m of the range of global IDs.
  int numGlobalIds() const {
    return mappedGlobalIds.size();
  }

  // Returns the size n of the range of local IDs.
  int numLocalIds() const {
    assert(!numMappedGlobalIds.empty());
    return numMappedGlobalIds.back();
  }

  // Returns true if globalId is mapped to a local ID.
  bool isGlobalIdMapped(const int globalId) const {
    return mappedGlobalIds[globalId];
  }

  // Returns the local ID globalId is mapped to.
  int toLocalId(const int globalId) const {
    assert(isGlobalIdMapped(globalId));
    return numMappedGlobalIdsBefore(globalId);
  }

  // Returns the number of mapped global IDs before globalId.
  int numMappedGlobalIdsBefore(const int globalId) const {
    assert(globalId >= 0); assert(globalId < numGlobalIds());
    int sum = numMappedGlobalIds[globalId / k];
    int blockIndex = globalId / k * k / BitVector::BITS_PER_BLOCK;
    for (; blockIndex < globalId / BitVector::BITS_PER_BLOCK; ++blockIndex)
      sum += bitCount(mappedGlobalIds.block(blockIndex));
    const int bitIndexWithinBlock = globalId % BitVector::BITS_PER_BLOCK;
    return sum + bitCountBeforeIndex(mappedGlobalIds.block(blockIndex), bitIndexWithinBlock);
  }

 private:
  const BitVector& mappedGlobalIds;        // The set of mapped global IDs.
  std::vector<int32_t> numMappedGlobalIds; // Stores the number of mapped IDs before every k-th ID.
};
