#pragma once

#include <array>
#include <cassert>
#include <type_traits>

#include <vectorclass/vectorclass.h>

#include "DataStructures/Labels/ParentInfo.h"
#include "Tools/Simd/HorizontalOps.h"
#include "Tools/TemplateProgramming.h"

// A set of consistent distance and parent labels for Dijkstra's algorithm. The template arguments
// specify the number of shortest paths computed simultaneously and the kind of parent information
// that should be collected.
template <int logSearches, ParentInfo parentInfo>
struct SimdLabelSet {
  static_assert(logSearches >= 2, "The SIMD label set requires at least 4 simultaneous searches.");

 public:
  // The number of simultaneous shortest-path computations.
  static constexpr int logK = logSearches;
  static constexpr int K = 1 << logK;

  // Vectors of multiple data elements for use with SIMD instructions.
  using BooleanVector = std::conditional_t<logK == 2, Vec4ib, Vec8ib>;
  using IntegerVector = std::conditional_t<logK == 2, Vec4i, Vec8i>;

  // Flags indicating whether parent vertices and/or parent edges should be collected.
  static constexpr bool KEEP_PARENT_VERTICES = parentInfo != ParentInfo::NO_PARENT_INFO;
  static constexpr bool KEEP_PARENT_EDGES    = parentInfo == ParentInfo::FULL_PARENT_INFO;

 private:
  static constexpr int VECTOR_SIZE = logK == 2 ? 4 : 8; // The number of data elements per vector.
  static constexpr int NUM_VECTORS = K / VECTOR_SIZE; // The number of vectors per label.

  // Arrays of vectors that store all k values of a packed label.
  using BooleanLabel = std::array<BooleanVector, NUM_VECTORS>;
  using IntegerLabel = std::array<IntegerVector, NUM_VECTORS>;

 public:
  // A mask that marks a subset of components in a packed distance label. For example, the result
  // of a less-than comparison between two multiple-source distance labels a and b is a mask that
  // indicates for which components i it holds that a[i] < b[i].
  class LabelMask {
   public:
    // Constructs an uninitialized mask.
    LabelMask() = default;

    // Constructs a mask with all k components set to val. Converting constructor.
    LabelMask(const bool val) {
      for (int i = 0; i < NUM_VECTORS; ++i)
        isMarked[i] = val;
    }

    // Takes the logical AND of this and the specified mask.
    LabelMask& operator&=(const LabelMask& rhs) {
      for (int i = 0; i < NUM_VECTORS; ++i)
        isMarked[i] &= rhs.isMarked[i];
      return *this;
    }

    // Returns a const reference to the i-th block of flags in this mask.
    const BooleanVector& operator[](const int i) const {
      assert(i >= 0); assert(i < NUM_VECTORS);
      return isMarked[i];
    }

    // Returns a reference to the i-th block of flags in this mask.
    BooleanVector& operator[](const int i) {
      assert(i >= 0); assert(i < NUM_VECTORS);
      return isMarked[i];
    }

    // Returns true if this mask marks at least one component.
    operator bool() const {
      BooleanVector tmp = isMarked[0];
      for (int i = 1; i < NUM_VECTORS; ++i)
        tmp |= isMarked[i];
      return horizontal_or(tmp);
    }

    BooleanLabel isMarked; // Flags indicating for each component if it is marked.
  };

  // A packed distance label for a vertex, storing k distance values. Each value maintains the
  // tentative distance from a different simultaneous source.
  class DistanceLabel {
   public:
    // A class simulating the behavior of references to a single distance value in a packed label.
    class Reference {
      // Only DistanceLabel is allowed to construct references to a single distance value.
      friend class DistanceLabel;

     public:
      // Sets the distance value to val.
      Reference& operator=(const int val) {
        block.insert(localIndex, val);
        return *this;
      }

      // Returns the distance value.
      operator int() const {
        return block.extract(localIndex);
      }

     private:
      // Constructs a reference to the i-th distance value in the specified block.
      Reference(IntegerVector& block, const int i) : block(block), localIndex(i) {}

      IntegerVector& block; // The block containing the distance value.
      const int localIndex; // The index of the distance value in the block above.
    };

    // Constructs an uninitialized distance label.
    DistanceLabel() = default;

    // Constructs a distance label with all k values set to val. Converting constructor.
    DistanceLabel(const int val) {
      for (int i = 0; i < NUM_VECTORS; ++i)
        values[i] = val;
    }

    // Returns a reference to the i-th distance value in this label.
    Reference operator[](const int i) {
      assert(i >= 0); assert(i < K);
      return Reference(values[i / VECTOR_SIZE], i % VECTOR_SIZE);
    }

    // Returns the packed sum of lhs and rhs.
    friend DistanceLabel operator+(const DistanceLabel& lhs, const DistanceLabel& rhs) {
      DistanceLabel sum;
      for (int i = 0; i < NUM_VECTORS; ++i)
        sum.values[i] = lhs.values[i] + rhs.values[i];
      return sum;
    }

    // Returns a mask that indicates for which components i it holds that lhs[i] < rhs[i].
    friend LabelMask operator<(const DistanceLabel& lhs, const DistanceLabel& rhs) {
      LabelMask mask;
      for (int i = 0; i < NUM_VECTORS; ++i)
        mask[i] = lhs.values[i] < rhs.values[i];
      return mask;
    }

    // Returns the priority of this label.
    int getKey() const {
      IntegerVector packedMin = values[0];
      for (int i = 1; i < NUM_VECTORS; ++i)
        packedMin = ::min(packedMin, values[i]);
      return horizontal_min(packedMin);
    }

    // Returns the horizontal maximum amongst the packed distance values in this label.
    int horizontalMax() const {
      IntegerVector packedMax = values[0];
      for (int i = 1; i < NUM_VECTORS; ++i)
        packedMax = ::max(packedMax, values[i]);
      return horizontal_max(packedMax);
    }

    // Take the packed minimum of this and the specified label.
    void min(const DistanceLabel& other) {
      for (int i = 0; i < NUM_VECTORS; ++i)
        values[i] = ::min(values[i], other.values[i]);
    }

   private:
    IntegerLabel values; // The k distance values, one for each simultaneous source.
  };

 private:
  // A packed label for a vertex, storing k parent edges.
  class ParentEdge {
   public:
    // Returns the parent edge on the shortest path from the i-th source.
    int edge(const int i) const {
      assert(i >= 0); assert(i < K);
      return edges[i / VECTOR_SIZE].extract(i % VECTOR_SIZE);
    }

    // Sets the parent edge to e on all shortest paths specified by mask.
    void setEdge(const int e, const LabelMask& mask) {
      for (int i = 0; i < NUM_VECTORS; ++i)
        edges[i] = select(mask[i], e, edges[i]);
    }

   private:
    IntegerLabel edges; // The k parent edges, one for each simultaneous source.
  };

 public:
  // A packed label for a vertex, storing k parent vertices and possibly k parent edges.
  class ParentLabel : public std::conditional_t<KEEP_PARENT_EDGES, ParentEdge, EmptyClass> {
   public:
    // Returns the parent vertex on the shortest path from the i-th source.
    int vertex(const int i) const {
      assert(i >= 0); assert(i < K);
      return vertices[i / VECTOR_SIZE].extract(i % VECTOR_SIZE);
    }

    // Sets the parent vertex to u on all shortest paths specified by mask.
    void setVertex(const int u, const LabelMask& mask) {
      for (int i = 0; i < NUM_VECTORS; ++i)
        vertices[i] = select(mask[i], u, vertices[i]);
    }

   private:
    IntegerLabel vertices; // The k parent vertices, one for each simultaneous source.
  };
};
