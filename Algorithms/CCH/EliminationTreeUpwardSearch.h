#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <vector>

#include "Algorithms/CH/CH.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Labels/Containers/ParentLabelContainer.h"
#include "DataStructures/Labels/Containers/SimpleDistanceLabelContainer.h"
#include "DataStructures/Queues/TournamentTree.h"
#include "Tools/Constants.h"

// An upward search in an elimination tree query. It enumerates all the vertices on a path in the
// elimination tree from a source vertex to the root, relaxing their outgoing edges. Depending on
// the label set, the algorithm keeps parent vertices and/or edges and computes multiple shortest
// paths simultaneously, possibly using SSE or AVX instructions.
template <typename LabelSetT>
class EliminationTreeUpwardSearch {
 private:
  using LabelMask = typename LabelSetT::LabelMask;         // Marks subset of components in label.
  using DistanceLabel = typename LabelSetT::DistanceLabel; // The distance label of a vertex.

 public:
  // Constructs an upward search instance.
  EliminationTreeUpwardSearch(
      const CH::SearchGraph& graph, const std::vector<int32_t>& eliminationTree,
      const DistanceLabel& tentativeDist)
      : searchGraph(graph),
        eliminationTree(eliminationTree),
        tentativeDistances(tentativeDist),
        distanceLabels(graph.numVertices()),
        parent(graph),
        nextVertices({INVALID_VERTEX}) {}

  // Initializes the labels of the source vertices.
  void init(const std::array<int, LabelSetT::K>& sources) {
    nextVertices.build(sources);
    distanceLabels[searchGraph.numVertices() - 1] = INFTY;
    for (int i = 0; i < LabelSetT::K; ++i) {
      const int s = sources[i];
      distanceLabels[s][i] = 0;
      parent.setVertex(s, INVALID_VERTEX, true);
      parent.setEdge(s, INVALID_EDGE, true);
    }
  }

  // Relaxes the edges out of the next vertex.
  void settleNextVertex() {
    const int u = getNextVertex();
    assert(u >= 0); assert(u < searchGraph.numVertices());
    DistanceLabel& distToU = distanceLabels[u];
#ifdef NO_FAST_ELIMINATION_TREE_QUERY
    if (distToU < INFTY)
#else
    if (distToU < tentativeDistances)
#endif
      FORALL_INCIDENT_EDGES(searchGraph, u, e) {
        const int v = searchGraph.edgeHead(e);
        DistanceLabel& distToV = distanceLabels[v];
        const DistanceLabel tentativeDist = distToU + searchGraph.traversalCost(e);
        const LabelMask mask = tentativeDist < distToV;
        if (mask) {
          distToV.min(tentativeDist);
          parent.setVertex(v, u, mask);
          parent.setEdge(v, e, mask);
        }
      }
    distToU = INFTY;

    // Find the next vertex to be settled. If two or more searches merged at u, block all but one.
    nextVertices.deleteMin(eliminationTree[u]);
    while (getNextVertex() == u)
      nextVertices.deleteMin(INFTY);
  }

  // Returns the distance label of the specified vertex.
  const DistanceLabel& getDistanceLabel(const int v) {
    return distanceLabels[v];
  }

  // Returns the vertices on the shortest path from the i-th source to t in reverse order.
  std::vector<int> getReversePath(int t, const int i) {
    return parent.getReversePath(t, i);
  }

  // Returns the edges on the shortest path from the i-th source to t in reverse order.
  std::vector<int> getReverseEdgePath(int t, const int i) {
    return parent.getReverseEdgePath(t, i);
  }

  // Returns the vertex to be settled next.
  int getNextVertex() const {
    return nextVertices.minKey();
  }

 private:
  using DistanceLabelCont = SimpleDistanceLabelContainer<DistanceLabel>;
  using ParentLabelCont = ParentLabelContainer<CH::SearchGraph, LabelSetT>;

  const CH::SearchGraph& searchGraph;      // The upward or downward graph.
  const std::vector<int>& eliminationTree; // eliminationTree[v] is the parent of v in the tree.
  const DistanceLabel& tentativeDistances; // One tentative distance for each simultaneous search.

  DistanceLabelCont distanceLabels;             // The distance labels of the vertices.
  ParentLabelCont parent;                       // The parent information for each vertex.
  TournamentTree<LabelSetT::logK> nextVertices; // The vertices settled next in the k searches.
};
