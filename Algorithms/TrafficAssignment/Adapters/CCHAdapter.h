#pragma once

#include <cassert>
#include <vector>

#include <routingkit/customizable_contraction_hierarchy.h>
#include <routingkit/nested_dissection.h>

#include "Algorithms/CH/CHConversion.h"
#include "Algorithms/CH/CHQuery.h"
#include "Algorithms/CH/ContractionHierarchy.h"
#include "DataStructures/Graph/Attributes/EdgeIdAttribute.h"
#include "DataStructures/Graph/Attributes/LatLngAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "DataStructures/Utilities/OriginDestination.h"

namespace trafficassignment {

// An adapter that makes CCHs usable in the all-or-nothing assignment procedure.
template <typename InputGraphT, typename WeightT>
class CCHAdapter {
 private:
  // The type of the CH resulting from perfectly customizing the CCH.
  using CHGraph = StaticGraph<VertexAttrs<>, EdgeAttrs<EdgeIdAttribute, WeightT>>;
  using CH = ContractionHierarchy<CHGraph, WeightT>;

 public:
  // The input graph type.
  using InputGraph = InputGraphT;

  // Constructs an adapter for CCHs.
  CCHAdapter(const InputGraph& graph) : inputGraph(graph), chSearch(perfectCH) {
    assert(graph.numEdges() > 0); assert(graph.isDefrag());
  }

  // Invoked before the first iteration.
  void preprocess() {
    // Convert the input graph to RoutingKit's graph representation.
    std::vector<float> lats(inputGraph.numVertices());
    std::vector<float> lngs(inputGraph.numVertices());
    std::vector<unsigned int> tails(inputGraph.numEdges());
    std::vector<unsigned int> heads(inputGraph.numEdges());
    FORALL_VERTICES(inputGraph, u) {
      lats[u] = inputGraph.latLng(u).latInDeg();
      lngs[u] = inputGraph.latLng(u).lngInDeg();
      FORALL_INCIDENT_EDGES(inputGraph, u, e) {
        tails[e] = u;
        heads[e] = inputGraph.edgeHead(e);
      }
    }

    // Compute a nested dissection order for the input graph.
    std::vector<unsigned int> order;
    order = RoutingKit::compute_nested_node_dissection_order_using_inertial_flow(
        inputGraph.numVertices(), tails, heads, lats, lngs);

    // Build the metric-independent CCH.
    cch = RoutingKit::CustomizableContractionHierarchy(order, tails, heads);
    const int* const weights = &inputGraph.template get<WeightT>(0);
    currentMetric = {cch, reinterpret_cast<const unsigned int*>(weights)};
  }

  // Invoked before each iteration.
  void customize() {
    // Customize the CCH using perfect witness searches.
    const int numEdges = inputGraph.numEdges();
    perfectCH = convert<CH>(
        currentMetric.build_contraction_hierarchy_using_perfect_witness_search(), numEdges);
    chSearch.resize();
  }

  // Computes the shortest path between the specified OD-pair.
  void query(const OriginDestination& od) {
    chSearch.run(perfectCH.rank(od.origin), perfectCH.rank(od.destination));
  }

  // Returns the length of the shortest path computed last.
  int getDistance(const int /*dst*/) {
    return chSearch.getDistance();
  }

  // Returns the edges on the (packed) shortest path computed last.
  std::vector<int> getPackedEdgePath(const int /*dst*/) {
    return chSearch.getPackedEdgePath();
  }

  // Returns the first constituent edge of shortcut s.
  int getShortcutsFirstEdge(const int s) const {
    return perfectCH.shortcutsFirstEdge(s);
  }

  // Returns the second constituent edge of shortcut s.
  int getShortcutsSecondEdge(const int s) const {
    return perfectCH.shortcutsSecondEdge(s);
  }

  // Returns the number of shortcut edges.
  int getNumShortcuts() const {
    return perfectCH.numShortcuts();
  }

 private:
  using CCH = RoutingKit::CustomizableContractionHierarchy;
  using CCHMetric = RoutingKit::CustomizableContractionHierarchyMetric;
  using CHQuery = StandardCHQuery<CH, BasicLabelSet<1, ParentInfo::FULL_PARENT_INFO>, false>;

  const InputGraph& inputGraph;     // The input graph.
  CCH cch;                          // The metric-independent CCH.
  CCHMetric currentMetric;          // The current metric for the CCH.
  CH perfectCH;                     // The CH resulting from perfectly customizing the CCH.
  CHQuery chSearch;                 // A CH search on the perfect CH.
};

}
