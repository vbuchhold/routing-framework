#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <utility>
#include <vector>

#include "Algorithms/CH/CH.h"
#include "Algorithms/CCH/UpwardEliminationTreeSearch.h"
#include "Algorithms/Dijkstra/Dijkstra.h"
#include "Tools/Constants.h"

// Implementation of an elimination tree query, which computes shortest paths in CCH without using
// priority queues. Depending on the used label set, it keeps parent vertices and/or edges, and
// computes multiple shortest paths simultaneously, optionally using SSE or AVX instructions.
template <typename LabelSetT>
class EliminationTreeQuery {
 private:
  using DistanceLabel = typename LabelSetT::DistanceLabel;
  using ParentLabel = typename LabelSetT::ParentLabel;

  // The pruning criterion for an elimination tree query that computes k shortest paths
  // simultaneously. We can prune the search at v if d_i(v) >= mu_i for all i = 1, ..., k.
  struct PruningCriterion {
    // Constructs a pruning criterion for an elimination tree query.
    PruningCriterion(const DistanceLabel& tentativeDistances) noexcept
        : tentativeDistances(&tentativeDistances) {}

    // Returns true if the search can be pruned at v.
    template <typename DistanceLabelContT>
    bool operator()(const int, const DistanceLabel& distToV, const DistanceLabelContT&) const {
      return !(distToV < *tentativeDistances);
    }

    const DistanceLabel* tentativeDistances; // One tentative distance per simultaneous search.
  };

  static constexpr int K = LabelSetT::K; // The number of simultaneous shortest-path computations.

 public:
  // Constructs an elimination tree query instance.
  EliminationTreeQuery(const CH& ch, const std::vector<int32_t>& eliminTree)
#ifdef NO_FAST_ELIMINATION_TREE_QUERY
      : forwardSearch(ch.upwardGraph(), eliminTree),
        reverseSearch(ch.downwardGraph(), eliminTree) {
#else
      : forwardSearch(ch.upwardGraph(), eliminTree, {tentativeDistances}),
        reverseSearch(ch.downwardGraph(), eliminTree, {tentativeDistances}) {
#endif
    assert(ch.upwardGraph().numVertices() == eliminTree.size());
  }

  // Move constructor.
  EliminationTreeQuery(EliminationTreeQuery&& other) noexcept
      : forwardSearch(std::move(other.forwardSearch)),
        reverseSearch(std::move(other.reverseSearch)) {
#ifndef NO_FAST_ELIMINATION_TREE_QUERY
    forwardSearch.pruneSearch = {tentativeDistances};
    reverseSearch.pruneSearch = {tentativeDistances};
#endif
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
    reverseSearch.distanceLabels[reverseSearch.searchGraph.numVertices() - 1] = INFTY;
    forwardSearch.init(sources);
    reverseSearch.init(targets);
    tentativeDistances = INFTY;
    while (forwardSearch.nextVertices.minKey() != INVALID_VERTEX)
      if (forwardSearch.nextVertices.minKey() <= reverseSearch.nextVertices.minKey()) {
        updateTentativeDistances(forwardSearch.nextVertices.minKey());
        forwardSearch.distanceLabels[forwardSearch.settleNextVertex()] = INFTY;
      } else {
        reverseSearch.distanceLabels[reverseSearch.settleNextVertex()] = INFTY;
      }
  }

  // Returns the length of the i-th shortest path.
  int getDistance(const int i = 0) {
    return tentativeDistances[i];
  }

  // Returns the edges in the upward graph on the up segment of the up-down path (in reverse order).
  const std::vector<int32_t>& getUpEdgePath(const int i = 0) {
    assert(tentativeDistances[i] != INFTY);
    return forwardSearch.getReverseEdgePath(meetingVertices.vertex(i), i);
  }

  // Returns the edges in the downward graph on the down segment of the up-down path.
  const std::vector<int32_t>& getDownEdgePath(const int i = 0) {
    assert(tentativeDistances[i] != INFTY);
    return reverseSearch.getReverseEdgePath(meetingVertices.vertex(i), i);
  }

 private:
  // Checks whether the path via v improves the tentative distance for any search.
  void updateTentativeDistances(const int v) {
    const auto distances = forwardSearch.distanceLabels[v] + reverseSearch.distanceLabels[v];
    meetingVertices.setVertex(v, distances < tentativeDistances);
    tentativeDistances.min(distances);
  }

  using UpwardSearch =
#ifdef NO_FAST_ELIMINATION_TREE_QUERY
      UpwardEliminationTreeSearch<LabelSetT>;
#else
      UpwardEliminationTreeSearch<LabelSetT, PruningCriterion>;
#endif

  UpwardSearch forwardSearch;       // The forward search from the source(s).
  UpwardSearch reverseSearch;       // The reverse search from the target(s).

  DistanceLabel tentativeDistances; // One tentative distance per simultaneous search.
  ParentLabel meetingVertices;      // One meeting vertex per simultaneous search.
};
