#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <limits>
#include <utility>
#include <vector>

#include "DataStructures/Graph/Graph.h"
#include "Tools/BinaryIO.h"
#include "Tools/Bitwise.h"
#include "Tools/Constants.h"

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

  // Returns a partition with the specified maximum cell size.
  std::vector<int> readOffPartition(const int maxCellSize) const {
    assert(maxCellSize > 0);
    const int numVertices = packedSideIds.size();
    std::vector<int> partition(numVertices, INVALID_ID);
    std::vector<std::pair<int, uint64_t>> ord(numVertices);
    for (int v = 0; v < numVertices; ++v)
      ord[v] = {v, packedSideIds[v]};
    std::sort(ord.begin(), ord.end(), [](const auto& u, const auto& v) {
      return u.second < v.second;
    });

    // Assign vertices to cells.
    int firstVertex = 0;
    int nextCellId = 0;
    for (int lastVertex = maxCellSize; lastVertex < numVertices; lastVertex += maxCellSize) {
      while (highestDifferingBit(ord[firstVertex].second, ord[lastVertex - 1].second) >
          highestDifferingBit(ord[lastVertex - 1].second, ord[lastVertex].second))
        --lastVertex;
      for (int i = firstVertex; i < lastVertex; ++i)
        partition[ord[i].first] = nextCellId;
      firstVertex = lastVertex;
      ++nextCellId;
    }
    for (int i = firstVertex; i < numVertices; ++i)
      partition[ord[i].first] = nextCellId;
    for (int i = 0; i < numVertices; ++i)
      assert(partition[i] != INVALID_ID);
    return partition;
  }

  // Returns a contraction order for the specified graph.
  template <typename GraphT>
  std::vector<int> readOffContractionOrder(const GraphT& graph) const {
    assert(graph.numVertices() == packedSideIds.size());
    // Compute for each vertex the level on which it is contained in a separator.
    std::vector<uint64_t> mask(graph.numVertices());
    const auto reverseGraph = graph.getReverseGraph();
    FORALL_VERTICES(graph, u) {
      uint64_t isSeparatorVertexOnLevel = 0;
      FORALL_INCIDENT_EDGES(graph, u, e)
        isSeparatorVertexOnLevel |= ~packedSideIds[u] & packedSideIds[graph.edgeHead(e)];
      FORALL_INCIDENT_EDGES(reverseGraph, u, e)
        isSeparatorVertexOnLevel |= ~packedSideIds[u] & packedSideIds[reverseGraph.edgeHead(e)];
      const int msb = highestOneBit(isSeparatorVertexOnLevel);
      mask[u] = msb == 63 ? 0 : uint64_t(-1) << (msb + 1);
    }

    std::vector<int> ord(graph.numVertices());
    FORALL_VERTICES(graph, u)
      ord[u] = u;

    // Compute the contraction order.
    std::sort(ord.begin(), ord.end(), [&mask, this](const auto u, const auto v) {
      const auto minMask = std::min(mask[u], mask[v]);
      return (packedSideIds[u] & minMask) < (packedSideIds[v] & minMask) ||
          ((packedSideIds[u] & minMask) == (packedSideIds[v] & minMask) &&
          (mask[u] > mask[v] || (mask[u] == mask[v] && u < v)));
    });
    return ord;
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
