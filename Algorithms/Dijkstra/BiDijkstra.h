#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <type_traits>
#include <vector>

#include "Tools/Constants.h"

namespace bidijkstra {

// The stopping criterion for a standard bidirectional Dijkstra search computing k shortest paths.
// We can stop the search as soon as mu_i <= Qf.minKey + Qr.minKey for all i = 1, ..., k.
template <typename QueueT>
struct BiDijkstraStoppingCriterion {
  // Constructs a stopping criterion for a standard bidirectional search.
  BiDijkstraStoppingCriterion(const QueueT& forwardQueue, const QueueT& reverseQueue,
                              const int& maxTentativeDistance)
      : forwardQueue(forwardQueue),
        reverseQueue(reverseQueue),
        maxTentativeDistance(maxTentativeDistance) {}

  // Returns true if we can stop the forward search.
  bool stopForwardSearch() const {
    return forwardQueue.empty() || reverseQueue.empty() ||
        maxTentativeDistance <= forwardQueue.minKey() + reverseQueue.minKey();
  }

  // Returns true if we can stop the reverse search.
  bool stopReverseSearch() const {
    return forwardQueue.empty() || reverseQueue.empty() ||
        maxTentativeDistance <= forwardQueue.minKey() + reverseQueue.minKey();
  }

  const QueueT& forwardQueue;      // The priority queue of the forward search.
  const QueueT& reverseQueue;      // The priority queue of the reverse search.
  const int& maxTentativeDistance; // The largest of all k tentative distances.
};

}

// Implementation of a bidirectional search. Depending on the underlying Dijkstra implementation, it
// keeps track of parent vertices and/or edges, and computes multiple shortest paths simultaneously,
// possibly using SSE or AVX instructions. The algorithm can be used with various stopping criteria,
// allowing it to be used as CH query algorithm.
template <typename DijkstraT,
          template <typename> class StoppingCriterionT = bidijkstra::BiDijkstraStoppingCriterion>
class BiDijkstra {
 private:
  using Graph = typename DijkstraT::Graph; // The graph type on which we compute shortest paths.

  static constexpr int K = DijkstraT::K; // The number of simultaneous shortest-path computations.

 public:
  // Constructs a bidirectional search instance.
  BiDijkstra(const Graph& graph, const Graph& reverseGraph,
             const typename DijkstraT::PruningCriterion& pruneForwardSearch = {},
             const typename DijkstraT::PruningCriterion& pruneReverseSearch = {})
      : forwardSearch(graph, pruneForwardSearch),
        reverseSearch(reverseGraph, pruneReverseSearch),
        stoppingCriterion(forwardSearch.queue, reverseSearch.queue, maxTentativeDistance) {}

  // Ensures that the internal data structures fit for the size of the graph.
  void resize() {
    forwardSearch.resize();
    reverseSearch.resize();
  }

  // Run a bidirectional search from s to t.
  void run(const int s, const int t) {
    std::array<int, K> sources;
    std::array<int, K> targets;
    std::fill(sources.begin(), sources.end(), s);
    std::fill(targets.begin(), targets.end(), t);
    run(sources, targets);
  }

  // Runs a bidirectional search that computes multiple shortest paths simultaneously.
  void run(const std::array<int, K>& sources, const std::array<int, K>& targets) {
    forwardSearch.init(sources);
    reverseSearch.init(targets);
    tentativeDistances = INFTY;
    maxTentativeDistance = INFTY;
    bool advanceForward = false;
    while (!stoppingCriterion.stopForwardSearch() || !stoppingCriterion.stopReverseSearch()) {
      advanceForward = !advanceForward; // Alternate between forward and reverse search.
      if ((advanceForward && !stoppingCriterion.stopForwardSearch()) ||
          stoppingCriterion.stopReverseSearch())
        updateTentativeDistances(forwardSearch.settleNextVertex());
      else
        updateTentativeDistances(reverseSearch.settleNextVertex());
    }
  }

  // Returns the length of the i-th shortest path.
  int getDistance(const int i = 0) {
    return tentativeDistances[i];
  }

  // Returns the vertices on the i-th shortest path.
  std::vector<int> getPath(const int i = 0) {
    assert(tentativeDistances[i] != INFTY);
    auto subpath1 = forwardSearch.getReversePath(meetingVertices.vertex(i), i);
    auto subpath2 = reverseSearch.getReversePath(meetingVertices.vertex(i), i);
    std::reverse(subpath1.begin(), subpath1.end());
    subpath1.pop_back();
    subpath1.insert(subpath1.end(), subpath2.begin(), subpath2.end());
    return subpath1;
  }

  // Returns the edges in the forward graph on the path to the meeting vertex (in reverse order).
  std::vector<int> getEdgePathToMeetingVertex(const int i = 0) {
    assert(tentativeDistances[i] != INFTY);
    return forwardSearch.getReverseEdgePath(meetingVertices.vertex(i), i);
  }

  // Returns the edges in the reverse graph on the path from the meeting vertex.
  std::vector<int> getEdgePathFromMeetingVertex(const int i = 0) {
    assert(tentativeDistances[i] != INFTY);
    return reverseSearch.getReverseEdgePath(meetingVertices.vertex(i), i);
  }

 private:
  using DistanceLabel = typename DijkstraT::DistanceLabel; // The distance label of a vertex.
  using ParentLabel = typename DijkstraT::ParentLabel;     // The parent label of a vertex.

  // Checks whether the path s-u-t improves the tentative distance for any search.
  void updateTentativeDistances(const int u) {
    DistanceLabel dist = forwardSearch.distanceLabels[u] + reverseSearch.distanceLabels[u];
    meetingVertices.setVertex(u, dist < tentativeDistances);
    tentativeDistances.min(dist);
    maxTentativeDistance = tentativeDistances.horizontalMax();
  }

  using StoppingCriterion = StoppingCriterionT<typename DijkstraT::Queue>;

  DijkstraT forwardSearch;             // The forward search from the source(s).
  DijkstraT reverseSearch;             // The reverse search from the target(s).
  StoppingCriterion stoppingCriterion; // The stopping criterion for the bidirectional search.

  DistanceLabel tentativeDistances; // One tentative distance for each simultaneous search.
  ParentLabel meetingVertices;      // One meeting vertex for each simultaneous search.
  int maxTentativeDistance;         // The largest of all k tentative distances.
};
