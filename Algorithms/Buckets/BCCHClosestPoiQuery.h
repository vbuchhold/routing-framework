#pragma once

#include "Algorithms/Buckets/BCCHBucketBuilder.h"
#include "Algorithms/Buckets/BucketClosestPoiQuery.h"
#include "Algorithms/CCH/CCH.h"
#include "Algorithms/CCH/UpwardEliminationTreeSearch.h"
#include "Algorithms/CH/CH.h"
#include "Algorithms/Dijkstra/Dijkstra.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"

// This class implements a BCCH-based closest-POI query in a road network. It uses an elimination
// tree search (which does not apply stall-on-demand) for the forward search.
class BCCHClosestPoiQuery : public BucketClosestPoiQuery<
    BCCHBucketBuilder,
    UpwardEliminationTreeSearch<
        BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>,
        dij::CompoundCriterion<
            buckets::StopWhenFrontierFartherThanClosestPois,
            buckets::ScanBucket>>> {
 public:
  // Creates a BCCH-based closest-POI query in the specified CCH.
  explicit BCCHClosestPoiQuery(const CCH& cch, const CH& minCH)
      : BucketClosestPoiQuery(
            minCH.upwardGraph().numVertices(),
            new BucketBuilder(minCH.upwardGraph(), minCH.downwardGraph(), cch.getEliminationTree()),
            ForwardCHSearch(
                minCH.upwardGraph(), cch.getEliminationTree(),
                {{currentlyClosestPois},
                    {numPoisToBeReported, bucketContainer, currentlyClosestPois}})) {}
};
