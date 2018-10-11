#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <vector>

#include "Algorithms/CH/CH.h"
#include "Algorithms/CCH/EliminationTreeUpwardSearch.h"
#include "Tools/Constants.h"

// An implementation of an elimination tree query, computing shortest paths in CCHs without using
// priority queues. Depending on the label set, the algorithm keeps parent vertices and/or edges,
// and computes multiple shortest paths simultaneously, possibly using SSE or AVX instructions.
template <typename LabelSetT>
class EliminationTreeQuery {
 public:
  static constexpr int K = LabelSetT::K; // The number of simultaneous shortest-path computations.

  // Constructs an elimination tree query instance.
  EliminationTreeQuery(const CH& ch, const std::vector<int32_t>& eliminationTree)
      : ch(ch),
        forwardSearch(ch.upwardGraph(), eliminationTree, tentativeDistances),
        reverseSearch(ch.downwardGraph(), eliminationTree, tentativeDistances) {}

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

  // Returns the vertices on the i-th up-down path.
  std::vector<int> getUpDownPath(const int i = 0) {
    assert(tentativeDistances[i] != INFTY);
    auto subpath1 = forwardSearch.getReversePath(meetingVertices.vertex(i), i);
    auto subpath2 = reverseSearch.getReversePath(meetingVertices.vertex(i), i);
    std::reverse(subpath1.begin(), subpath1.end());
    subpath1.pop_back();
    subpath1.insert(subpath1.end(), subpath2.begin(), subpath2.end());
    return subpath1;
  }

  // Returns the edges in the upward graph on the up part of the up-down path (in reverse order).
  std::vector<int> getUpEdgePath(const int i = 0) {
    assert(tentativeDistances[i] != INFTY);
    return forwardSearch.getReverseEdgePath(meetingVertices.vertex(i), i);
  }

  // Returns the edges in the downward graph on the down part of the up-down path.
  std::vector<int> getDownEdgePath(const int i = 0) {
    assert(tentativeDistances[i] != INFTY);
    return reverseSearch.getReverseEdgePath(meetingVertices.vertex(i), i);
  }

  // Returns the edges in the input graph on the i-th shortest path.
  std::vector<int> getEdgePath(const int i = 0) {
    auto upPath = getUpEdgePath(i);
    auto downPath = getDownEdgePath(i);
    std::reverse(downPath.begin(), downPath.end());
    std::for_each(downPath.begin(), downPath.end(), [](int& e) { e = -e - 1; });
    downPath.insert(downPath.end(), upPath.begin(), upPath.end());

    std::vector<int> fullPath;
    const auto& upGraph = ch.upwardGraph();
    const auto& downGraph = ch.downwardGraph();

    while (!downPath.empty()) {
      const auto e = downPath.back();
      downPath.pop_back();
      const auto& unpackInfo = e >= 0 ? upGraph.unpackingInfo(e) : downGraph.unpackingInfo(-e - 1);
      if (unpackInfo.second == INVALID_EDGE) {
        fullPath.push_back(unpackInfo.first);
      } else {
        downPath.push_back(unpackInfo.second);
        downPath.push_back(-unpackInfo.first - 1);
      }
    }
    return fullPath;
  }

 private:
  using DistanceLabel = typename LabelSetT::DistanceLabel; // The distance label of a vertex.
  using ParentLabel = typename LabelSetT::ParentLabel;     // The parent information for a vertex.

  // Checks if the path s-u-t improves the tentative distance for any search.
  void updateTentativeDistances(const int u) {
    DistanceLabel dist = forwardSearch.getDistanceLabel(u) + reverseSearch.getDistanceLabel(u);
    meetingVertices.setVertex(u, dist < tentativeDistances);
    tentativeDistances.min(dist);
  }

  using UpwardSearch = EliminationTreeUpwardSearch<LabelSetT>;

  const CH& ch;                     // The CH on which we compute shortest paths.
  UpwardSearch forwardSearch;       // The forward search from the source vertices.
  UpwardSearch reverseSearch;       // The reverse search from the target vertices.
  DistanceLabel tentativeDistances; // One tentative distance for each simultaneous search.
  ParentLabel meetingVertices;      // One meeting vertex for each simultaneous search.
};
