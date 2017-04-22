#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <type_traits>

#include "DataStructures/Labels/ParentInfo.h"
#include "Tools/TemplateProgramming.h"

// A set of consistent distance and parent labels for Dijkstra's algorithm. The template arguments
// specify the number of shortest paths computed simultaneously and the kind of parent information
// that should be collected.
template <int logSearches, ParentInfo parentInfo>
struct BasicLabelSet {
 public:
  // The number of simultaneous shortest-path computations.
  static constexpr int logK = logSearches;
  static constexpr int K = 1 << logK;

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

    // Constructs a mask with all k components set to val. Converting constructor.
    LabelMask(const bool val) {
      std::fill(isMarked.begin(), isMarked.end(), val);
    }

    // Takes the logical AND of this and the specified mask.
    LabelMask& operator&=(const LabelMask& rhs) {
      for (int i = 0; i < K; ++i)
        isMarked[i] &= rhs.isMarked[i];
      return *this;
    }

    // Returns the i-th block of flags in this mask.
    bool operator[](const int i) const {
      assert(i >= 0); assert(i < K);
      return isMarked[i];
    }

    // Returns a reference to the i-th block of flags in this mask.
    bool& operator[](const int i) {
      assert(i >= 0); assert(i < K);
      return isMarked[i];
    }

    // Returns true if this mask marks at least one component.
    operator bool() const {
      bool res = isMarked[0];
      for (int i = 1; i < K; ++i)
        res |= isMarked[i];
      return res;
    }

    std::array<bool, K> isMarked; // Flags indicating for each component if it is marked.
  };

  // A packed distance label for a vertex, storing k distance values. Each value maintains the
  // tentative distance from a different simultaneous source.
  class DistanceLabel {
   public:
    // Constructs an uninitialized distance label.
    DistanceLabel() = default;

    // Constructs a distance label with all k values set to val. Converting constructor.
    DistanceLabel(const int val) {
      std::fill(values.begin(), values.end(), val);
    }

    // Returns a reference to the i-th distance value in this label.
    int& operator[](const int i) {
      assert(i >= 0); assert(i < K);
      return values[i];
    }

    // Returns the packed sum of lhs and rhs.
    friend DistanceLabel operator+(const DistanceLabel& lhs, const DistanceLabel& rhs) {
      DistanceLabel sum;
      for (int i = 0; i < K; ++i)
        sum.values[i] = lhs.values[i] + rhs.values[i];
      return sum;
    }

    // Returns a mask that indicates for which components i it holds that lhs[i] < rhs[i].
    friend LabelMask operator<(const DistanceLabel& lhs, const DistanceLabel& rhs) {
      LabelMask mask;
      for (int i = 0; i < K; ++i)
        mask[i] = lhs.values[i] < rhs.values[i];
      return mask;
    }

    // Returns the priority of this label.
    int getKey() const {
      int min = values[0];
      for (int i = 1; i < K; ++i)
        min = std::min(min, values[i]);
      return min;
    }

    // Returns the horizontal maximum amongst the packed distance values in this label.
    int horizontalMax() const {
      int max = values[0];
      for (int i = 1; i < K; ++i)
        max = std::max(max, values[i]);
      return max;
    }

    // Take the packed minimum of this and the specified label.
    void min(const DistanceLabel& other) {
      for (int i = 0; i < K; ++i)
        values[i] = std::min(values[i], other.values[i]);
    }

   private:
    std::array<int, K> values; // The k distance values, one for each simultaneous source.
  };

 private:
  // A packed label for a vertex, storing k parent edges.
  class ParentEdge {
   public:
    // Returns the parent edge on the shortest path from the i-th source.
    int edge(const int i) const {
      assert(i >= 0); assert(i < K);
      return edges[i];
    }

    // Sets the parent edge to e on all shortest paths specified by mask.
    void setEdge(const int e, const LabelMask& mask) {
      for (int i = 0; i < K; ++i)
        edges[i] = mask[i] ? e : edges[i];
    }

   private:
    std::array<int, K> edges; // The k parent edges, one for each simultaneous source.
  };

 public:
  // A packed label for a vertex, storing k parent vertices and possibly k parent edges.
  class ParentLabel : public std::conditional_t<KEEP_PARENT_EDGES, ParentEdge, EmptyClass> {
   public:
    // Returns the parent vertex on the shortest path from the i-th source.
    int vertex(const int i) const {
      assert(i >= 0); assert(i < K);
      return vertices[i];
    }

    // Sets the parent vertex to u on all shortest paths specified by mask.
    void setVertex(const int u, const LabelMask& mask) {
      for (int i = 0; i < K; ++i)
        vertices[i] = mask[i] ? u : vertices[i];
    }

   private:
    std::array<int, K> vertices; // The k parent vertices, one for each simultaneous source.
  };
};
