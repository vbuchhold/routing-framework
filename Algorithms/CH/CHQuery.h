#pragma once

#include <algorithm>
#include <array>
#include <type_traits>
#include <vector>

#include "Algorithms/CH/CH.h"
#include "Algorithms/Dijkstra/BiDijkstra.h"
#include "Algorithms/Dijkstra/Dijkstra.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Labels/Containers/StampedDistanceLabelContainer.h"
#include "DataStructures/Queues/AddressableKHeap.h"
#include "Tools/Constants.h"

// Implementation of a CH query. Depending on the used label set, it keeps track of parent vertices
// and/or edges, and computes multiple shortest paths simultaneously, possibly using SSE or AVX
// instructions. The algorithm can be used with different distance label containers and queues.
template <
    template <typename> class DistanceLabelContT, typename LabelSetT, typename QueueT,
    bool useStallOnDemand = true>
class CHQuery {
 public:
  static constexpr int K = LabelSetT::K; // The number of simultaneous shortest-path computations.

  // Constructs a CH search instance that uses stall-on-demand.
  template <bool cond = useStallOnDemand>
  CHQuery(std::enable_if_t<cond, const CH&> ch)
      : ch(ch),
        chSearch(
            ch.upwardGraph(), ch.downwardGraph(),
            CHQueryPruningCriterion(ch.downwardGraph()),
            CHQueryPruningCriterion(ch.upwardGraph())) {}

  // Constructs a CH search instance that does not use stall-on-demand.
  template <bool cond = useStallOnDemand>
  CHQuery(std::enable_if_t<!cond, const CH&> ch)
      : ch(ch), chSearch(ch.upwardGraph(), ch.downwardGraph()) {}

  // Runs a CH search from s to t.
  void run(const int s, const int t) {
    chSearch.run(s, t);
  }

  // Runs a CH search that computes multiple shortest paths simultaneously.
  void run(const std::array<int, K>& sources, const std::array<int, K>& targets) {
    chSearch.run(sources, targets);
  }

  // Returns the length of the i-th shortest path.
  int getDistance(const int i = 0) {
    return chSearch.getDistance(i);
  }

  // Returns the vertices on the i-th up-down path.
  std::vector<int> getUpDownPath(const int i = 0) {
    return chSearch.getPath(i);
  }

  // Returns the edges in the upward graph on the up part of the up-down path (in reverse order).
  std::vector<int> getUpEdgePath(const int i = 0) {
    return chSearch.getEdgePathToMeetingVertex(i);
  }

  // Returns the edges in the downward graph on the down part of the up-down path.
  std::vector<int> getDownEdgePath(const int i = 0) {
    return chSearch.getEdgePathFromMeetingVertex(i);
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
  // The pruning criterion for a CH query, also known as stall-on-demand.
  struct CHQueryPruningCriterion {
    using LabelMask = typename LabelSetT::LabelMask;     // The label mask type we use.
    using DistLabel = typename LabelSetT::DistanceLabel; // The distance label type we use.

    // Constructs a pruning criterion for a CH query.
    CHQueryPruningCriterion(const CH::SearchGraph& oppositeGraph) : oppositeGraph(oppositeGraph) {}

    // Returns true if the search can be pruned at u.
    template <typename DistLabelContainerT>
    bool operator()(const int u, const DistLabel& distToU, DistLabelContainerT& distLabels) const {
      LabelMask continueSearch = true;
      FORALL_INCIDENT_EDGES(oppositeGraph, u, e)
        continueSearch &=
            distToU < distLabels[oppositeGraph.edgeHead(e)] + oppositeGraph.traversalCost(e);
      return !continueSearch;
    }

    const CH::SearchGraph& oppositeGraph; // The down graph if we prune the up search or vice versa.
  };

  // The stopping criterion for a CH query computing k shortest paths simultaneously. We can stop
  // the forward search as soon as mu_i <= Qf.minKey for all i = 1, ..., k. The reverse search is
  // stopped in the same way.
  template <typename>
  struct CHQueryStoppingCriterion {
    // Constructs a stopping criterion for a CH query.
    CHQueryStoppingCriterion(const QueueT& forwardQueue, const QueueT& reverseQueue,
                             const int& maxTentativeDistance)
        : forwardQueue(forwardQueue),
          reverseQueue(reverseQueue),
          maxTentativeDistance(maxTentativeDistance) {}

    // Returns true if we can stop the forward search.
    bool stopForwardSearch() const {
      return forwardQueue.empty() || maxTentativeDistance <= forwardQueue.minKey();
    }

    // Returns true if we can stop the reverse search.
    bool stopReverseSearch() const {
      return reverseQueue.empty() || maxTentativeDistance <= reverseQueue.minKey();
    }

    const QueueT& forwardQueue;      // The priority queue of the forward search.
    const QueueT& reverseQueue;      // The priority queue of the reverse search.
    const int& maxTentativeDistance; // The largest of all k tentative distances.
  };

  using CHDijkstraStall = Dijkstra<
      CH::SearchGraph, TraversalCostAttribute, DistanceLabelContT, LabelSetT, QueueT,
      CHQueryPruningCriterion>;
  using CHDijkstraNoStall = Dijkstra<
      CH::SearchGraph, TraversalCostAttribute, DistanceLabelContT, LabelSetT, QueueT>;
  using CHDijkstra = std::conditional_t<useStallOnDemand, CHDijkstraStall, CHDijkstraNoStall>;

  const CH& ch;                                              // The CH on which we work.
  BiDijkstra<CHDijkstra, CHQueryStoppingCriterion> chSearch; // The modified bidirectional search.
};

// An alias template for a standard CH search.
template <typename LabelSetT, bool useStallOnDemand = true>
using StandardCHQuery =
    CHQuery<StampedDistanceLabelContainer, LabelSetT, AddressableQuadheap, useStallOnDemand>;
