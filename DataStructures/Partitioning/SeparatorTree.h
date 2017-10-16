#pragma once

#include <cassert>
#include <cstdint>
#include <fstream>
#include <limits>
#include <vector>

#include "Tools/BinaryIO.h"
#include "Tools/Bitwise.h"

// A separator tree of a graph computed by recursive bisection. It stores a bit for each vertex and
// level that indicates the side which the vertex belongs to. All bits for a vertex are packed in a
// single 64-bit integer with the higher bits representing the higher levels. Separator trees allow
// efficiently reading off (multilevel) partitions and contraction orders.
class SeparatorTree {
 public:
  // Constructs a separator tree with all side IDs initialized to zero.
  explicit SeparatorTree(const int numVertices) : packedSideIds(numVertices) {}

  // Constructs a separator tree from a binary file.
  explicit SeparatorTree(std::ifstream& in) {
    readFrom(in);
  }

  // Returns the side that contains v on the specified level. The levels are indexed top-down.
  bool getSide(const int v, const int level) const {
    assert(v >= 0); assert(v < packedSideIds.size());
    return getBit(packedSideIds[v], std::numeric_limits<uint64_t>::digits - level - 1);
  }

  // Sets the side that contains v on the specified level to val. The levels are indexed top-down.
  void setSide(const int v, const int level, const bool val) {
    assert(v >= 0); assert(v < packedSideIds.size());
    setBit(packedSideIds[v], std::numeric_limits<uint64_t>::digits - level - 1, val);
  }

  // Read a separator tree from a binary file.
  void readFrom(std::ifstream& in) {
    bio::read(in, packedSideIds);
  }

  // Write a separator tree to a binary file.
  void writeTo(std::ofstream& out) const {
    bio::write(out, packedSideIds);
  }

 private:
  std::vector<uint64_t> packedSideIds; // The packed side IDs for each vertex.
};
