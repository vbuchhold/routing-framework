#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <utility>
#include <vector>

#include "Tools/Constants.h"

namespace bidij {

// The stopping criterion for a standard bidirectional search that computes k shortest paths. We can
// stop the search as soon as mu_i <= Qf.minKey + Qr.minKey for all i = 1, ..., k.
template <typename QueueT>
struct StoppingCriterion {
  // Constructs a stopping criterion for a standard bidirectional search.
  StoppingCriterion(
      const QueueT& forwardQueue, const QueueT& reverseQueue, const int& maxTentativeDistance)
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
// keeps parent vertices and/or edges, and computes multiple shortest paths simultaneously,
// optionally using SSE or AVX instructions. The algorithm can be used with different stopping
// criteria, allowing it to be used as CH query algorithm.
template <
    typename DijkstraT, template <typename> class StoppingCriterionT = bidij::StoppingCriterion>
class BiDijkstra {
 private:
  static constexpr int K = DijkstraT::K; // The number of simultaneous shortest-path computations.

 public:
  using Graph = typename DijkstraT::Graph; // The graph on which we compute shortest paths.

  // Constructs a bidirectional search instance.
  BiDijkstra(
      const Graph& forwardGraph, const Graph& reverseGraph,
      typename DijkstraT::PruningCriterion pruneForwardSearch = {},
      typename DijkstraT::PruningCriterion pruneReverseSearch = {})
      : forwardSearch(forwardGraph, {}, pruneForwardSearch),
        reverseSearch(reverseGraph, {}, pruneReverseSearch),
        stoppingCriterion(forwardSearch.queue, reverseSearch.queue, maxTentativeDistance) {
    assert(forwardGraph.numVertices() == reverseGraph.numVertices());
  }

  // Move constructor.
  BiDijkstra(BiDijkstra&& other) noexcept
      : forwardSearch(std::move(other.forwardSearch)),
        reverseSearch(std::move(other.reverseSearch)),
        stoppingCriterion(forwardSearch.queue, reverseSearch.queue, maxTentativeDistance),
        tentativeDistances(other.tentativeDistances),
        meetingVertices(other.meetingVertices),
        maxTentativeDistance(other.maxTentativeDistance) {}

  // Runs a bidirectional search from s to t.
  void run(const int s, const int t) {
    std::array<int, K> sources;
    std::array<int, K> targets;
    sources.fill(s);
    targets.fill(t);
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
      advanceForward = !advanceForward; // Alternate between the forward and reverse search.
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

  // Returns the edges in the forward graph on the path to the meeting vertex (in reverse order).
  const std::vector<int32_t>& getEdgePathToMeetingVertex(const int i = 0) {
    assert(tentativeDistances[i] != INFTY);
    return forwardSearch.getReverseEdgePath(meetingVertices.vertex(i), i);
  }

  // Returns the edges in the reverse graph on the path from the meeting vertex.
  const std::vector<int32_t>& getEdgePathFromMeetingVertex(const int i = 0) {
    assert(tentativeDistances[i] != INFTY);
    return reverseSearch.getReverseEdgePath(meetingVertices.vertex(i), i);
  }

 private:
  // Checks whether the path via v improves the tentative distance for any search.
  void updateTentativeDistances(const int v) {
    const auto distances = forwardSearch.distanceLabels[v] + reverseSearch.distanceLabels[v];
    meetingVertices.setVertex(v, distances < tentativeDistances);
    tentativeDistances.min(distances);
    maxTentativeDistance = tentativeDistances.horizontalMax();
  }

  using StoppingCriterion = StoppingCriterionT<typename DijkstraT::Queue>;
  using DistanceLabel = typename DijkstraT::DistanceLabel;
  using ParentLabel = typename DijkstraT::ParentLabel;

  DijkstraT forwardSearch;             // The forward search from the source(s).
  DijkstraT reverseSearch;             // The reverse search from the target(s).
  StoppingCriterion stoppingCriterion; // The criterion used to stop the search.

  DistanceLabel tentativeDistances; // One tentative distance per simultaneous search.
  ParentLabel meetingVertices;      // One meeting vertex per simultaneous search.
  int maxTentativeDistance;         // The largest of all k tentative distances.
};
