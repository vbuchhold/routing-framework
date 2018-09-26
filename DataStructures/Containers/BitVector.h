#pragma once

#include <cassert>
#include <limits>
#include <vector>

#include "Tools/Simd/AlignedVector.h"
#include "Tools/Bitwise.h"

// This class implements a vector of bits. The bits are divided into 64-bit blocks. Different bits
// can be modified concurrently when they are in different blocks. The first block is aligned on a
// cache line boundary.
class BitVector {
 public:
  using Block = uint64_t; // The unsigned integer type in which the bits are stored.

  // A class that simulates the behavior of references to a single bit in a bit vector.
  class Reference {
    friend BitVector; // Only bit vectors are allowed to construct a reference.

   public:
    // Returns the value of the bit.
    operator bool() const {
      return getBit(block, bitIndex);
    }

    // Sets the bit to bitValue.
    Reference& operator=(const bool bitValue) {
      setBit(block, bitIndex, bitValue);
      return *this;
    }

    // Sets the bit to the value of the specified bit.
    Reference& operator=(const Reference& other) {
      return *this = static_cast<bool>(other);
    }

   private:
    // Constructs a reference to the bit with the specified index in the specified block.
    Reference(Block& block, const int bitIndex) : block(block), bitIndex(bitIndex) {
      assert(bitIndex >= 0); assert(bitIndex < std::numeric_limits<Block>::digits);
    }

    Block& block;       // The block in which the bit is stored.
    const int bitIndex; // The index of the bit within its block.
  };

  // The number of bits in a block.
  static constexpr int BITS_PER_BLOCK = std::numeric_limits<Block>::digits;

  // Constructs a bit vector of the specified size. All bits are initialized to init.
  explicit BitVector(const int size, const bool init = false)
      : blocks((size + BITS_PER_BLOCK - 1) / BITS_PER_BLOCK, init ? -1 : 0), numBits(size) {}

  // Returns the number of bits in this bit vector.
  int size() const {
    return numBits;
  }

  // Returns the number of blocks in this bit vector.
  int numBlocks() const {
    return blocks.size();
  }

  // Returns the bit with the specified index.
  bool operator[](const int bitIndex) const {
    assert(bitIndex >= 0); assert(bitIndex < size());
    return getBit(blocks[bitIndex / BITS_PER_BLOCK], bitIndex % BITS_PER_BLOCK);
  }

  // Returns a reference to the bit with the specified index.
  Reference operator[](const int bitIndex) {
    assert(bitIndex >= 0); assert(bitIndex < size());
    return {blocks[bitIndex / BITS_PER_BLOCK], bitIndex % BITS_PER_BLOCK};
  }

  // Returns the block with the specified index.
  Block block(const int blockIndex) const {
    assert(blockIndex >= 0); assert(blockIndex < numBlocks());
    return blocks[blockIndex];
  }

  // Returns the index of the first one-bit. If no such bit exists then -1 is returned.
  int firstSetBit() const {
    int blockIndex = 0;
    while (blocks[blockIndex] == 0 && blockIndex < blocks.size()) ++blockIndex;
    if (blockIndex == blocks.size()) return -1;
    return blockIndex * BITS_PER_BLOCK + numTrailingZeros(blocks[blockIndex]);
  }

  // Returns the index of the first one-bit that occurs after fromIndex.
  // If no such bit exists then -1 is returned.
  int nextSetBit(int fromIndex) const {
    assert(fromIndex >= 0); assert(fromIndex < size());
    if (++fromIndex == size()) return -1;
    int blockIndex = fromIndex / BITS_PER_BLOCK;
    const auto firstBlock = blocks[blockIndex] >> fromIndex % BITS_PER_BLOCK;
    if (firstBlock != 0) return fromIndex + numTrailingZeros(firstBlock);
    ++blockIndex;
    while (blocks[blockIndex] == 0 && blockIndex < blocks.size()) ++blockIndex;
    if (blockIndex == blocks.size()) return -1;
    return blockIndex * BITS_PER_BLOCK + numTrailingZeros(blocks[blockIndex]);
  }

 private:
  AlignedVector<Block> blocks; // The blocks in which the bits are stored.
  int numBits;                 // The number of bits in this bit vector.
};
