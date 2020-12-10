#pragma once

#include <cassert>

#include "Algorithms/Buckets/BucketBuilder.h"
#include "Algorithms/CH/CH.h"
#include "Algorithms/CH/CHQuery.h"
#include "Algorithms/Dijkstra/Dijkstra.h"
#include "DataStructures/Graph/Attributes/TraversalCostAttribute.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"

// This class implements the construction of static standard-CH buckets for a given set of targets.
// It uses a Dijkstra-based CH search that applies the stall-on-demand optimization.
class BCHBucketBuilder : public BucketBuilder<Dijkstra<
    CH::SearchGraph, TraversalCostAttribute, BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>,
    dij::NoCriterion,
    dij::CompoundCriterion<
        CHQuery<BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>>::PruningCriterion,
        buckets::CollectSearchSpace>>> {
 public:
  // Creates a facility for building standard-CH buckets.
  BCHBucketBuilder(const CH::SearchGraph& upwardGraph, const CH::SearchGraph& downwardGraph)
      : BucketBuilder(
            upwardGraph.numVertices(),
            CHSearch(downwardGraph, {}, {{upwardGraph}, {targetId, searchSpaces}})) {
    assert(upwardGraph.numVertices() == downwardGraph.numVertices());
  }
};
