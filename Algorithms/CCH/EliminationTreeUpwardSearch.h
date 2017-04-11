#pragma once

#include <cassert>
#include <vector>

#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Labels/Containers/ParentLabelContainer.h"
#include "DataStructures/Labels/Containers/SimpleDistanceLabelContainer.h"
#include "Tools/Constants.h"

// An upward search in an elimination tree query. It enumerates all the vertices on a path in the
// elimination tree from a source vertex to the root, relaxing their outgoing edges. Depending on
// the label set, the algorithm keeps parent vertices and/or edges.
template <typename SearchGraphT, typename WeightT, typename LabelSetT>
class EliminationTreeUpwardSearch {
 private:
  using LabelMask = typename LabelSetT::LabelMask;         // Marks subset of components in label.
  using DistanceLabel = typename LabelSetT::DistanceLabel; // The distance label of a vertex.
  using ParentLabel   = typename LabelSetT::ParentLabel;   // The parent information for a vertex.

 public:
  // Constructs an upward search instance.
  EliminationTreeUpwardSearch(
      const SearchGraphT& graph, const std::vector<int>& eliminationTree,
      const DistanceLabel& tentativeDist)
      : searchGraph(graph),
        eliminationTree(eliminationTree),
        tentativeDistances(tentativeDist),
        distanceLabels(graph.numVertices()),
        parent(graph),
        nextVertex(INVALID_VERTEX) {}

  // Ensures that the internal data structures fit the size of the graph.
  void resize() {
    distanceLabels.resize(searchGraph.numVertices());
    parent.resize();
  }

  // Initializes the labels of the source vertex.
  void init(const int s) {
    distanceLabels[s][0] = 0;
    parent.setVertex(s, INVALID_VERTEX, LabelMask(0));
    parent.setEdge(s, INVALID_EDGE, LabelMask(0));
    nextVertex = s;
  }

  // Relaxes the edges out of the next vertex.
  void settleNextVertex() {
    assert(nextVertex >= 0); assert(nextVertex < searchGraph.numVertices());
    DistanceLabel& distToU = distanceLabels[nextVertex];
    if (distToU < tentativeDistances)
      FORALL_INCIDENT_EDGES(searchGraph, nextVertex, e) {
        const int v = searchGraph.edgeHead(e);
        DistanceLabel& distToV = distanceLabels[v];
        const DistanceLabel tentativeDist = distToU + searchGraph.template get<WeightT>(e);
        const LabelMask mask = tentativeDist < distToV;
        if (mask) {
          distToV.min(tentativeDist);
          parent.setVertex(v, nextVertex, mask);
          parent.setEdge(v, e, mask);
        }
      }
    distToU = INFTY;
    nextVertex = eliminationTree[nextVertex];
  }

  // Returns the distance label of the specified vertex.
  const DistanceLabel& getDistanceLabel(const int v) {
    return distanceLabels[v];
  }

  // Returns the vertices on the shortest path to t in reverse order.
  std::vector<int> getReversePath(int t) {
    return parent.getReversePath(t);
  }

  // Returns the edges on the shortest path to t in reverse order.
  std::vector<int> getReverseEdgePath(int t) {
    return parent.getReverseEdgePath(t);
  }

  // Returns the vertex to be settled next.
  int getNextVertex() const {
    return nextVertex;
  }

 private:
  using DistanceLabelCont = SimpleDistanceLabelContainer<DistanceLabel>;
  using ParentLabelCont = ParentLabelContainer<SearchGraphT, LabelSetT>;

  const SearchGraphT& searchGraph;         // The upward or downward graph.
  const std::vector<int>& eliminationTree; // eliminationTree[v] is the parent of v in the tree.
  const DistanceLabel& tentativeDistances; // One tentative distance for each simultaneous search.

  DistanceLabelCont distanceLabels; // The distance labels of the vertices.
  ParentLabelCont parent;           // The parent information for each vertex.
  int nextVertex;                   // The vertex to be settled next.
};
