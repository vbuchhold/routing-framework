#pragma once

#include <algorithm>
#include <cassert>
#include <vector>

#include "Algorithms/CCH/EliminationTreeUpwardSearch.h"
#include "Tools/Constants.h"

// An implementation of an elimination tree query, computing shortest paths in CCHs without using
// priority queues. Depending on the label set, the algorithm keeps parent vertices and/or edges.
template <typename CH, typename LabelSetT>
class EliminationTreeQuery {
 public:
  // Constructs an elimination tree query instance.
  EliminationTreeQuery(const CH& ch, const std::vector<int>& eliminationTree)
      : forwardSearch(ch.getUpwardGraph(), eliminationTree, tentativeDistances),
        reverseSearch(ch.getDownwardGraph(), eliminationTree, tentativeDistances) {}

  // Ensures that the internal data structures fit the size of the graph.
  void resize() {
    forwardSearch.resize();
    reverseSearch.resize();
  }

  // Runs an elimination tree query from s to t.
  void run(const int s, const int t) {
    forwardSearch.init(s);
    reverseSearch.init(t);
    tentativeDistances = INFTY;

    // Settles all vertices on the path in the elimination tree from s (t) to the LCA.
    while (forwardSearch.getNextVertex() != reverseSearch.getNextVertex())
      if (forwardSearch.getNextVertex() < reverseSearch.getNextVertex())
        forwardSearch.settleNextVertex();
      else
        reverseSearch.settleNextVertex();

    // Settles all vertices on the path in the elimination tree from the LCA to the root.
    while (forwardSearch.getNextVertex() != INVALID_VERTEX) {
      assert(forwardSearch.getNextVertex() == reverseSearch.getNextVertex());
      updateTentativeDistances(forwardSearch.getNextVertex());
      forwardSearch.settleNextVertex();
      reverseSearch.settleNextVertex();
    }
  }

  // Returns the shortest-path distance.
  int getDistance() {
    return tentativeDistances[0];
  }

  // Returns the vertices on the packed shortest path.
  std::vector<int> getPackedPath() {
    assert(tentativeDistances[0] != INFTY);
    std::vector<int> subpath1 = forwardSearch.getReversePath(meetingVertices.vertex(0));
    std::vector<int> subpath2 = reverseSearch.getReversePath(meetingVertices.vertex(0));
    std::reverse(subpath1.begin(), subpath1.end());
    subpath1.pop_back();
    subpath1.insert(subpath1.end(), subpath2.begin(), subpath2.end());
    return subpath1;
  }

  // Returns the edges on the packed shortest path.
  std::vector<int> getPackedEdgePath() {
    assert(tentativeDistances[0] != INFTY);
    std::vector<int> subpath1 = forwardSearch.getReverseEdgePath(meetingVertices.vertex(0));
    std::vector<int> subpath2 = reverseSearch.getReverseEdgePath(meetingVertices.vertex(0));
    std::reverse(subpath1.begin(), subpath1.end());
    subpath1.insert(subpath1.end(), subpath2.begin(), subpath2.end());
    return subpath1;
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

  UpwardSearch forwardSearch;       // The forward search from the source vertex.
  UpwardSearch reverseSearch;       // The reverse search from the target vertex.
  DistanceLabel tentativeDistances; // One tentative distance for each simultaneous search.
  ParentLabel meetingVertices;      // One meeting vertex for each simultaneous search.
};
