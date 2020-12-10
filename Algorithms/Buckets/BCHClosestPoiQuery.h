#pragma once

#include "Algorithms/Buckets/BCHBucketBuilder.h"
#include "Algorithms/Buckets/BucketClosestPoiQuery.h"
#include "Algorithms/CH/CH.h"
#include "Algorithms/CH/CHQuery.h"
#include "Algorithms/Dijkstra/Dijkstra.h"
#include "DataStructures/Graph/Attributes/TraversalCostAttribute.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"

// This class implements a BCH-based closest-POI query in a road network. It uses a Dijkstra-based
// CH search (which applies the stall-on-demand optimization) for the forward search.
class BCHClosestPoiQuery : public BucketClosestPoiQuery<
    BCHBucketBuilder,
    Dijkstra<
        CH::SearchGraph, TraversalCostAttribute, BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>,
        buckets::StopWhenFrontierFartherThanClosestPois,
        dij::CompoundCriterion<
            CHQuery<BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>>::PruningCriterion,
            buckets::ScanBucket>>> {
 public:
  // Creates a BCH-based closest-POI query in the specified CH.
  explicit BCHClosestPoiQuery(const CH& ch)
      : BucketClosestPoiQuery(
            ch.upwardGraph().numVertices(),
            new BucketBuilder(ch.upwardGraph(), ch.downwardGraph()),
            ForwardCHSearch(
                ch.upwardGraph(), {currentlyClosestPois},
                {{ch.downwardGraph()},
                    {numPoisToBeReported, bucketContainer, currentlyClosestPois}})) {}
};
