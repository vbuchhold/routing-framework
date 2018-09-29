#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>

#include "Algorithms/CCH/CCH.h"
#include "DataStructures/Graph/Graph.h"
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
    #pragma omp parallel
    #pragma omp single nowait
    computeCustomizedMetric();
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
  void runPerfectCustomization() noexcept {
    cch.forEachVertexTopDown([&](const int u) {
      FORALL_INCIDENT_EDGES(cch.getUpwardGraph(), u, lower) {
        const int v = cch.getUpwardGraph().edgeHead(lower);
        cch.forEachUpperTriangle(u, v, lower, [&](int, const int inter, const int upper) {
          if (upWeights[inter] + downWeights[upper] < upWeights[lower])
            upWeights[lower] = upWeights[inter] + downWeights[upper];
          if (upWeights[lower] + upWeights[upper] < upWeights[inter])
            upWeights[inter] = upWeights[lower] + upWeights[upper];
          if (upWeights[upper] + downWeights[inter] < downWeights[lower])
            downWeights[lower] = upWeights[upper] + downWeights[inter];
          if (downWeights[upper] + downWeights[lower] < downWeights[inter])
            downWeights[inter] = downWeights[upper] + downWeights[lower];
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
