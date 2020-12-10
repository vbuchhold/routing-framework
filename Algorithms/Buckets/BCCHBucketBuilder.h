#pragma once

#include <cassert>
#include <cstdint>
#include <vector>

#include "Algorithms/Buckets/BucketBuilder.h"
#include "Algorithms/CCH/UpwardEliminationTreeSearch.h"
#include "Algorithms/CH/CH.h"
#include "Algorithms/CH/CHQuery.h"
#include "Algorithms/Dijkstra/Dijkstra.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"

// This class implements the construction of static customizable-CH buckets for a given set of
// targets. It uses an elimination tree search that applies the stall-on-demand optimization.
class BCCHBucketBuilder : public BucketBuilder<UpwardEliminationTreeSearch<
    BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>,
    dij::CompoundCriterion<
        elimintree::PruningCriterion,
        CHQuery<BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>>::PruningCriterion,
        buckets::CollectSearchSpace>>> {
 public:
  // Creates a facility for building customizable-CH buckets.
  BCCHBucketBuilder(
      const CH::SearchGraph& upwardGraph, const CH::SearchGraph& downwardGraph,
      const std::vector<int32_t>& elimTree)
      : BucketBuilder(
            upwardGraph.numVertices(),
            CHSearch(downwardGraph, elimTree, {{}, {upwardGraph}, {targetId, searchSpaces}})) {
    assert(upwardGraph.numVertices() == downwardGraph.numVertices());
    assert(upwardGraph.numVertices() == elimTree.size());
  }
};
