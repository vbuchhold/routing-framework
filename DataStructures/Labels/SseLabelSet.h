#pragma once

#include <array>
#include <cassert>
#include <type_traits>

#include <emmintrin.h>
#include <smmintrin.h>
#include <xmmintrin.h>

#include "DataStructures/Labels/ParentInfo.h"
#include "Tools/Simd/Sse.h"
#include "Tools/TemplateProgramming.h"

// A set of consistent distance and parent labels for Dijkstra's algorithm. The template arguments
// specify the number of shortest paths computed simultaneously and the kind of parent information
// that should be collected.
template <int numSources, ParentInfo parentInfo>
struct SseLabelSet {
 public:
  static constexpr int K = numSources; // The number of simultaneous shortest-path computations.

  // Flags indicating whether parent vertices and/or parent edges should be collected.
  static constexpr bool KEEP_PARENT_VERTICES = parentInfo != ParentInfo::NO_PARENT_INFO;
  static constexpr bool KEEP_PARENT_EDGES    = parentInfo == ParentInfo::FULL_PARENT_INFO;

  // A mask that marks a subset of components in a packed distance label. For example, the result
  // of a less-than comparison between two multiple-source distance labels a and b is a mask that
  // indicates for which components i it holds that a[i] < b[i].
  class LabelMask {
   public:
    // Constructs an uninitialized mask.
    LabelMask() = default;

    // Constructs a mask that marks only the i-th component.
    explicit LabelMask(const int i) {
      assert(i >= 0); assert(i < K);
      for (int j = 0; j < K / 4; ++j)
        isMarked[j] = _mm_setzero_si128();
      isMarked[i / 4] = mm_insert_epi32(isMarked[i / 4], -1, i % 4);
    }

    // Constructs a mask with all k components set to val. Converting constructor.
    LabelMask(const bool val) {
      for (int i = 0; i < K / 4; ++i)
        isMarked[i] = _mm_set1_epi32(val * -1);
    }

    // Takes the logical AND of this and the specified mask.
    LabelMask& operator&=(const LabelMask& rhs) {
      for (int i = 0; i < K / 4; ++i)
        isMarked[i] = _mm_and_si128(isMarked[i], rhs.isMarked[i]);
      return *this;
    }

    // Returns the i-th block of flags in this mask.
    __m128i operator[](const int i) const {
      assert(i >= 0); assert(i < K / 4);
      return isMarked[i];
    }

    // Returns a reference to the i-th block of flags in this mask.
    __m128i& operator[](const int i) {
      assert(i >= 0); assert(i < K / 4);
      return isMarked[i];
    }

    // Returns true if this mask marks at least one component.
    operator bool() const {
      bool res = _mm_movemask_epi8(isMarked[0]);
      for (int i = 0; i < K / 4; ++i)
        res |= _mm_movemask_epi8(isMarked[i]);
      return res;
    }

    std::array<__m128i, K / 4> isMarked; // Flags indicating for each component if it is marked.
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
        block = mm_insert_epi32(block, val, localIndex);
        return *this;
      }

      // Returns the distance value.
      operator int() const {
        return mm_extract_epi32(block, localIndex);
      }

     private:
      // Constructs a reference to the i-th distance value in the specified block.
      Reference(__m128i& block, const int i) : block(block), localIndex(i) {}

      __m128i& block;       // The block containing the distance value.
      const int localIndex; // The index of the distance value in the block above.
    };

    // Constructs an uninitialized distance label.
    DistanceLabel() = default;

    // Constructs a distance label with all k values set to val. Converting constructor.
    DistanceLabel(const int val) {
      for (int i = 0; i < K / 4; ++i)
        values[i] = _mm_set1_epi32(val);
    }

    // Returns a reference to the i-th distance value in this label.
    Reference operator[](const int i) {
      assert(i >= 0); assert(i < K);
      return Reference(values[i / 4], i % 4);
    }

    // Returns the packed sum of this label plus rhs.
    DistanceLabel operator+(const DistanceLabel& rhs) const {
      DistanceLabel sum;
      for (int i = 0; i < K / 4; ++i)
        sum.values[i] = _mm_add_epi32(values[i], rhs.values[i]);
      return sum;
    }

    // Returns a mask that indicates for which components i it holds that lhs[i] < rhs[i].
    friend LabelMask operator<(const DistanceLabel& lhs, const DistanceLabel& rhs) {
      LabelMask mask;
      for (int i = 0; i < K / 4; ++i)
        mask[i] = _mm_cmpgt_epi32(rhs.values[i], lhs.values[i]);
      return mask;
    }

    // Returns the priority of this label.
    int getKey() const {
      __m128i min = values[0];
      for (int i = 1; i < K / 4; ++i)
        min = _mm_min_epi32(min, values[i]);

      // Horizontally compute the minimum amongst the distance values in min.
      min = _mm_min_epi32(min, _mm_shuffle_epi32(min, _MM_SHUFFLE(1, 0, 3, 2)));
      min = _mm_min_epi32(min, _mm_shuffle_epi32(min, _MM_SHUFFLE(2, 3, 0, 1)));
      return _mm_extract_epi32(min, 0);
    }

    // Returns the horizontal maximum amongst the packed distance values in this label.
    int horizontalMax() const {
      __m128i max = values[0];
      for (int i = 1; i < K / 4; ++i)
        max = _mm_max_epi32(max, values[i]);

      // Horizontally compute the maximum amongst the distance values in max.
      max = _mm_max_epi32(max, _mm_shuffle_epi32(max, _MM_SHUFFLE(1, 0, 3, 2)));
      max = _mm_max_epi32(max, _mm_shuffle_epi32(max, _MM_SHUFFLE(2, 3, 0, 1)));
      return _mm_extract_epi32(max, 0);
    }

    // Take the packed minimum of this and the specified label.
    void min(const DistanceLabel& other) {
      for (int i = 0; i < K / 4; ++i)
        values[i] = _mm_min_epi32(values[i], other.values[i]);
    }

   private:
    std::array<__m128i, K / 4> values; // The k distance values, one for each simultaneous source.
  };

 private:
  // A packed label for a vertex, storing k parent edges.
  class ParentEdge {
   public:
    // Returns the parent edge on the shortest path from the i-th source.
    int edge(const int i) const {
      assert(i >= 0); assert(i < K);
      return mm_extract_epi32(edges[i / 4], i % 4);
    }

    // Sets the parent edge to e on all shortest paths specified by mask.
    void setEdge(const int e, const LabelMask& mask) {
      for (int i = 0; i < K / 4; ++i)
        edges[i] = _mm_blendv_epi8(edges[i], _mm_set1_epi32(e), mask[i]);
    }

   private:
    std::array<__m128i, K / 4> edges; // The k parent edges, one for each simultaneous source.
  };

 public:
  // A packed label for a vertex, storing k parent vertices and possibly k parent edges.
  class ParentLabel : public std::conditional_t<KEEP_PARENT_EDGES, ParentEdge, EmptyClass> {
   public:
    // Returns the parent vertex on the shortest path from the i-th source.
    int vertex(const int i) const {
      assert(i >= 0); assert(i < K);
      return mm_extract_epi32(vertices[i / 4], i % 4);
    }

    // Sets the parent vertex to u on all shortest paths specified by mask.
    void setVertex(const int u, const LabelMask& mask) {
      for (int i = 0; i < K / 4; ++i)
        vertices[i] = _mm_blendv_epi8(vertices[i], _mm_set1_epi32(u), mask[i]);
    }

   private:
    std::array<__m128i, K / 4> vertices; // The k parent vertices, one for each simultaneous source.
  };
};
