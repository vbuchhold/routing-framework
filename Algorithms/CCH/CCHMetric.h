#pragma once

#include <algorithm>
#include <cstdint>
#include <utility>
#include <vector>

#include "Algorithms/CCH/CCH.h"
#include "Algorithms/CH/CH.h"
#include "DataStructures/Containers/ConcurrentLocalIdMap.h"
#include "DataStructures/Graph/Attributes/TraversalCostAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "Tools/Simd/AlignedVector.h"
#include "Tools/ConcurrentHelpers.h"
#include "Tools/Constants.h"

// This class encodes the actual cost of the edges in a customizable contraction hierarchy. It
// stores the edge weights and contains several sequential and parallel customization algorithms.
class CCHMetric {
 public:
  // Constructs an individual metric incorporating the specified input weights in the specified CCH.
  CCHMetric(const CCH& cch, const int32_t* const inputWeights)
      : cch(cch), inputWeights(inputWeights) {
    assert(inputWeights != nullptr);
    upWeights.resize(cch.getUpwardGraph().numEdges());
    downWeights.resize(cch.getUpwardGraph().numEdges());
  }

  // Incorporates the current input weights in this metric.
  void customize() {
    computeRespectingMetric();
    computeCustomizedMetric();
  }

  // Runs the perfect customization algorithm.
  void runPerfectCustomization() noexcept {
    runPerfectCustomization([](const int /*e*/) {}, [](const int /*e*/) {});
  }

  // Returns a weighted CH having the smallest possible number of edges for the given order.
  CH buildMinimumWeightedCH() {
    std::vector<int8_t> keepUpEdge;
    std::vector<int8_t> keepDownEdge;

    #pragma omp parallel sections
    {
      #pragma omp section
      keepUpEdge.resize(upWeights.size() + 1, true);
      #pragma omp section
      keepDownEdge.resize(upWeights.size() + 1, true);
    }

    keepUpEdge.back() = false;
    keepDownEdge.back() = false;
    customize();
    runPerfectCustomization(
        [&](const int e) { keepUpEdge[e] = false; },
        [&](const int e) { keepDownEdge[e] = false; });

    const auto& cchGraph = cch.getUpwardGraph();
    ConcurrentLocalIdMap<4> upEdgeIdMap(keepUpEdge);
    ConcurrentLocalIdMap<4> downEdgeIdMap(keepDownEdge);
    const auto numUpEdges = upEdgeIdMap.numLocalIds();
    const auto numDownEdges = downEdgeIdMap.numLocalIds();

    AlignedVector<CH::SearchGraph::OutEdgeRange> upOutEdges;
    AlignedVector<CH::SearchGraph::OutEdgeRange> downOutEdges;
    AlignedVector<int32_t> upEdgeHeads;
    AlignedVector<int32_t> downEdgeHeads;
    AlignedVector<TraversalCostAttribute::Type> upEdgeWeights;
    AlignedVector<TraversalCostAttribute::Type> downEdgeWeights;

    #pragma omp parallel sections
    {
      #pragma omp section
      upOutEdges.resize(cchGraph.numVertices() + 1);
      #pragma omp section
      downOutEdges.resize(cchGraph.numVertices() + 1);
      #pragma omp section
      upEdgeHeads.resize(numUpEdges);
      #pragma omp section
      downEdgeHeads.resize(numDownEdges);
      #pragma omp section
      upEdgeWeights.resize(numUpEdges);
      #pragma omp section
      downEdgeWeights.resize(numDownEdges);
    }

    #pragma omp parallel for schedule(dynamic, 2048)
    FORALL_VERTICES(cchGraph, v) {
      upOutEdges[v].first() = upEdgeIdMap.numMappedGlobalIdsBefore(cchGraph.firstEdge(v));
      downOutEdges[v].first() = downEdgeIdMap.numMappedGlobalIdsBefore(cchGraph.firstEdge(v));
    }

    #pragma omp parallel for schedule(dynamic, 2048)
    FORALL_EDGES(cchGraph, e) {
      if (keepUpEdge[e]) {
        const auto newIdx = upEdgeIdMap.toLocalId(e);
        upEdgeHeads[newIdx] = cchGraph.edgeHead(e);
        upEdgeWeights[newIdx] = upWeights[e];
      }
      if (keepDownEdge[e]) {
        const auto newIdx = downEdgeIdMap.toLocalId(e);
        downEdgeHeads[newIdx] = cchGraph.edgeHead(e);
        downEdgeWeights[newIdx] = downWeights[e];
      }
    }

    upOutEdges.back().first() = numUpEdges;
    downOutEdges.back().first() = numDownEdges;

    CH::SearchGraph upGraph(
        std::move(upOutEdges), std::move(upEdgeHeads), numUpEdges,
        std::move(upEdgeWeights));
    CH::SearchGraph downGraph(
        std::move(downOutEdges), std::move(downEdgeHeads), numDownEdges,
        std::move(downEdgeWeights));

    Permutation order;
    Permutation ranks;

    #pragma omp parallel sections
    {
      #pragma omp section
      order = cch.getContractionOrder();
      #pragma omp section
      ranks = cch.getRanks();
    }

    return {std::move(upGraph), std::move(downGraph), std::move(order), std::move(ranks)};
  }

 private:
  // Computes a respecting metric.
  void computeRespectingMetric() {
    upWeights.resize(cch.getUpwardGraph().numEdges());
    downWeights.resize(cch.getUpwardGraph().numEdges());
    #pragma omp parallel for schedule(static)
    FORALL_EDGES(cch.getUpwardGraph(), e) {
      upWeights[e] = INFTY;
      downWeights[e] = INFTY;
      cch.forEachUpwardInputEdge(e, [&](const int inputEdge) {
        if (inputWeights[inputEdge] < upWeights[e])
          upWeights[e] = inputWeights[inputEdge];
      });
      cch.forEachDownwardInputEdge(e, [&](const int inputEdge) {
        if (inputWeights[inputEdge] < downWeights[e])
          downWeights[e] = inputWeights[inputEdge];
      });
    }
  }

  // Computes a customized metric given a respecting one.
  void computeCustomizedMetric() noexcept {
    #pragma omp parallel
    #pragma omp single nowait
    if (omp_get_num_threads() == 1)
      computeCustomizedMetricSequentially();
    else
      computeCustomizedMetricInParallel();
  }

  // Computes a customized metric sequentially.
  void computeCustomizedMetricSequentially() noexcept {
    cch.forEachVertexBottomUp([&](const int u) {
      FORALL_INCIDENT_EDGES(cch.getUpwardGraph(), u, lower) {
        const int v = cch.getUpwardGraph().edgeHead(lower);
        cch.forEachUpperTriangle(u, v, lower, [&](int, const int inter, const int upper) {
          if (downWeights[lower] + upWeights[inter] < upWeights[upper])
            upWeights[upper] = downWeights[lower] + upWeights[inter];
          if (downWeights[inter] + upWeights[lower] < downWeights[upper])
            downWeights[upper] = downWeights[inter] + upWeights[lower];
          return true;
        });
      }
    });
  }

  // Computes a customized metric in parallel.
  void computeCustomizedMetricInParallel() noexcept {
    cch.forEachVertexBottomUp([&](const int u) {
      FORALL_INCIDENT_EDGES(cch.getUpwardGraph(), u, lower) {
        const int v = cch.getUpwardGraph().edgeHead(lower);
        cch.forEachUpperTriangle(u, v, lower, [&](int, const int inter, const int upper) {
          atomicFetchMin(upWeights[upper], downWeights[lower] + upWeights[inter]);
          atomicFetchMin(downWeights[upper], downWeights[inter] + upWeights[lower]);
          return true;
        });
      }
    });
  }

  // Runs the perfect customization algorithm.
  template <typename T1, typename T2>
  void runPerfectCustomization(T1 markUpEdgeForRemoval, T2 markDownEdgeForRemoval) noexcept {
    #pragma omp parallel
    #pragma omp single nowait
    cch.forEachVertexTopDown([&](const int u) {
      FORALL_INCIDENT_EDGES(cch.getUpwardGraph(), u, lower) {
        const int v = cch.getUpwardGraph().edgeHead(lower);
        cch.forEachUpperTriangle(u, v, lower, [&](int, const int inter, const int upper) {
          if (upWeights[inter] + downWeights[upper] < upWeights[lower]) {
            upWeights[lower] = upWeights[inter] + downWeights[upper];
            markUpEdgeForRemoval(lower);
          }
          if (upWeights[lower] + upWeights[upper] < upWeights[inter]) {
            upWeights[inter] = upWeights[lower] + upWeights[upper];
            markUpEdgeForRemoval(inter);
          }
          if (upWeights[upper] + downWeights[inter] < downWeights[lower]) {
            downWeights[lower] = upWeights[upper] + downWeights[inter];
            markDownEdgeForRemoval(lower);
          }
          if (downWeights[upper] + downWeights[lower] < downWeights[inter]) {
            downWeights[inter] = downWeights[upper] + downWeights[lower];
            markDownEdgeForRemoval(inter);
          }
          return true;
        });
      }
    });
  }

  const CCH& cch;                    // The associated CCH.
  const int32_t* const inputWeights; // The weights of the input edges.

  std::vector<int32_t> upWeights;   // The upward weights of the edges in the CCH.
  std::vector<int32_t> downWeights; // The downward weights of the edges in the CCH.
};
