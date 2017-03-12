#pragma once

#include <cassert>
#include <vector>

#include <routingkit/contraction_hierarchy.h>

#include "Algorithms/CH/CHConversion.h"
#include "Algorithms/CH/ContractionHierarchy.h"
#include "DataStructures/Graph/Graph.h"

// An implementation of the CH preprocessing. The vertices in the upward and downward graph are
// reordered by rank to improve data locality during queries.
template <typename InputGraphT, typename SearchGraphT, template <typename> class GetWeightT>
class CHPreprocessing {
 public:
  using CH = ContractionHierarchy<SearchGraphT, GetWeightT>; // The CH type we build.

  // Constructs a CH preprocessing instance.
  CHPreprocessing(const InputGraphT& inputGraph) : inputGraph(inputGraph) {
    assert(inputGraph.isDefrag());
  }

  // Contracts the input graph.
  CH run() {
    // Convert the input graph into RoutingKit's graph representation.
    std::vector<unsigned int> tails(inputGraph.numEdges());
    std::vector<unsigned int> heads(inputGraph.numEdges());
    std::vector<unsigned int> weights(inputGraph.numEdges());
    FORALL_VALID_EDGES(inputGraph, u, e) {
      tails[e] = u;
      heads[e] = inputGraph.edgeHead(e);
      weights[e] = getWeight(inputGraph, e);
    }

    const RoutingKit::ContractionHierarchy ch =
        RoutingKit::ContractionHierarchy::build(inputGraph.numVertices(), tails, heads, weights);
    return convert<CH>(ch, inputGraph.numEdges());
  }

 private:
  const InputGraphT& inputGraph;     // The input graph that should be contracted.
  GetWeightT<InputGraphT> getWeight; // A functor returning the edge weight used for routing.
};
