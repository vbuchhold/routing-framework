#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <vector>

#include "Algorithms/CCH/EliminationTreeUpwardSearch.h"
#include "Tools/Constants.h"

// An implementation of an elimination tree query, computing shortest paths in CCHs without using
// priority queues. Depending on the label set, the algorithm keeps parent vertices and/or edges,
// and computes multiple shortest paths simultaneously, possibly using SSE or AVX instructions.
template <typename CH, typename LabelSetT>
class EliminationTreeQuery {
 public:
  static constexpr int K = LabelSetT::K; // The number of simultaneous shortest-path computations.

  // Constructs an elimination tree query instance.
  EliminationTreeQuery(const CH& ch, const std::vector<int>& eliminationTree)
      : ch(ch),
        forwardSearch(ch.getUpwardGraph(), eliminationTree, tentativeDistances),
        reverseSearch(ch.getDownwardGraph(), eliminationTree, tentativeDistances) {}

  // Ensures that the internal data structures fit the size of the graph.
  void resize() {
    forwardSearch.resize();
    reverseSearch.resize();
  }

  // Runs an elimination tree query from s to t.
  void run(const int s, const int t) {
    std::array<int, K> sources;
    std::array<int, K> targets;
    sources.fill(s);
    targets.fill(t);
    run(sources, targets);
  }

  // Runs an elimination tree query that computes multiple shortest paths simultaneously.
  void run(const std::array<int, K>& sources, const std::array<int, K>& targets) {
    forwardSearch.init(sources);
    reverseSearch.init(targets);
    tentativeDistances = INFTY;
    while (forwardSearch.getNextVertex() != INVALID_VERTEX)
      if (forwardSearch.getNextVertex() <= reverseSearch.getNextVertex()) {
        updateTentativeDistances(forwardSearch.getNextVertex());
        forwardSearch.settleNextVertex();
      } else {
        reverseSearch.settleNextVertex();
      }
  }

  // Returns the length of the i-th shortest path.
  int getDistance(const int i = 0) {
    return tentativeDistances[i];
  }

  // Returns the vertices on the i-th packed shortest path.
  std::vector<int> getPackedPath(const int i = 0) {
    assert(tentativeDistances[i] != INFTY);
    std::vector<int> subpath1 = forwardSearch.getReversePath(meetingVertices.vertex(i), i);
    std::vector<int> subpath2 = reverseSearch.getReversePath(meetingVertices.vertex(i), i);
    std::reverse(subpath1.begin(), subpath1.end());
    subpath1.pop_back();
    subpath1.insert(subpath1.end(), subpath2.begin(), subpath2.end());
    return subpath1;
  }

  // Returns the edges on the i-th packed shortest path.
  std::vector<int> getPackedEdgePath(const int i = 0) {
    assert(tentativeDistances[i] != INFTY);
    std::vector<int> subpath1 = forwardSearch.getReverseEdgePath(meetingVertices.vertex(i), i);
    std::vector<int> subpath2 = reverseSearch.getReverseEdgePath(meetingVertices.vertex(i), i);
    std::reverse(subpath1.begin(), subpath1.end());
    subpath1.insert(subpath1.end(), subpath2.begin(), subpath2.end());
    return subpath1;
  }

  // Returns the edges on the i-th shortest path.
  std::vector<int> getEdgePath(const int i = 0) {
    auto packedPath = getPackedEdgePath(i);
    std::reverse(packedPath.begin(), packedPath.end());
    std::vector<int> fullPath;
    while (!packedPath.empty()) {
      const int e = packedPath.back();
      packedPath.pop_back();
      if (ch.isEdgeShortcut(e)) {
        packedPath.push_back(ch.shortcutsSecondEdge(e));
        packedPath.push_back(ch.shortcutsFirstEdge(e));
      } else {
        fullPath.push_back(e);
      }
    }
    return fullPath;
  }

 private:
  using DistanceLabel = typename LabelSetT::DistanceLabel; // The distance label of a vertex.
  using ParentLabel = typename LabelSetT::ParentLabel;     // The parent information for a vertex.

  // Checks whether the path s-u-t improves the tentative distance for any search.
  void updateTentativeDistances(const int u) {
    DistanceLabel dist = forwardSearch.getDistanceLabel(u) + reverseSearch.getDistanceLabel(u);
    meetingVertices.setVertex(u, dist < tentativeDistances);
    tentativeDistances.min(dist);
  }

  using SearchGraph = typename CH::SearchGraph;
  using Weight      = typename CH::Weight;
  using UpwardSearch = EliminationTreeUpwardSearch<SearchGraph, Weight, LabelSetT>;

  const CH& ch;                     // The CH on which we compute shortest paths.
  UpwardSearch forwardSearch;       // The forward search from the source vertices.
  UpwardSearch reverseSearch;       // The reverse search from the target vertices.
  DistanceLabel tentativeDistances; // One tentative distance for each simultaneous search.
  ParentLabel meetingVertices;      // One meeting vertex for each simultaneous search.
};
