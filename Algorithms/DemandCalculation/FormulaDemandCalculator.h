#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

#include "Algorithms/Dijkstra/Dijkstra.h"
#include "DataStructures/Graph/Attributes/TravelTimeAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "Tools/CommandLine/ProgressBar.h"
#include "Tools/Constants.h"
#include "Tools/OpenMP.h"
#include "Tools/Timer.h"

// A travel demand calculator based on evaluating radiation model's formula for each pair of
// vertices. For each vertex s in the graph, we run a Dijkstra search without early stopping.
// Whenever we settle a vertex u, we know the intervening opportunities between s and u, and
// we evaluate the formula for s and u.
template <typename GraphT>
class FormulaDemandCalculator {
 public:
  // Constructs a travel demand calculator for the specified network.
  explicit FormulaDemandCalculator(const GraphT& graph, const int seed, const bool verbose) noexcept
      : graph(graph), totPop(0), totPoi(1), seed(seed), verbose(verbose) {
    FORALL_VERTICES(graph, v) {
      totPop += graph.population(v);
      totPoi += graph.numOpportunities(v);
    }
    assert(totPop > 0);
    assert(totPoi > 1);
    assert(seed >= 0);
  }

  // Generates OD pairs and writes them to the specified file.
  void calculateDemand(
      int numODPairs, double lambda, double swapProb, const std::string& fileName) const {
    Timer timer;
    const auto maxNumSources = std::min(graph.numVertices(), DC_MAX_NUM_SOURCES);
    if (verbose) std::cout << "Calculating demand: ";
    ProgressBar bar(maxNumSources);

    assert(lambda >= 0); assert(lambda <= 1);
    using LabelSet = BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>;
    using Dijkstra = StandardDijkstra<GraphT, TravelTimeAttribute, LabelSet>;

    #pragma omp parallel
    {
      Dijkstra dijkstra(graph);
      std::minstd_rand rand(seed + omp_get_thread_num() + 1);
      std::bernoulli_distribution swapDist(swapProb);

      std::ofstream out(fileName + ".part" + std::to_string(omp_get_thread_num()));
      assert(out.good());

      #pragma omp for schedule(static, 1) nowait
      for (auto src = 0; src < maxNumSources; ++src) {
        auto srcPop = graph.population(src);           // The source population.
        auto srcPoi = graph.numOpportunities(src) + 1; // The number of opportunities at the source.
        auto intPoi = 0;                               // The number of intervening opportunities.
        if (srcPop == 0)
          continue;

        auto outflow = std::binomial_distribution<>(numODPairs, srcPop / totPop)(rand);
        auto normConst = 1 / (1 -
            (1 - std::pow(lambda, totPoi)) / (1 - std::pow(lambda, srcPoi)) * (srcPoi / totPoi));
        auto n = (1 - std::pow(lambda, srcPoi)) / srcPoi;

        dijkstra.init({{src}});
        dijkstra.settleNextVertex();
        while (!dijkstra.queue.empty()) {
          auto dst = dijkstra.settleNextVertex();
          auto dstPoi = graph.numOpportunities(dst);
          if (dstPoi == 0)
            continue;

          // Evaluate radiation model's formula for src and dst.
          auto ns = (1 - std::pow(lambda, srcPoi + intPoi)) / (srcPoi + intPoi);
          auto nns = (1 - std::pow(lambda, srcPoi + dstPoi + intPoi)) / (srcPoi + dstPoi + intPoi);
          auto prob = (ns - nns) / n;
          auto numTravelers = std::binomial_distribution<>(outflow, normConst * prob)(rand);

          for (auto i = 0; i < numTravelers; ++i)
            if (swapDist(rand))
              out << dst << ',' << src << '\n';
            else
              out << src << ',' << dst << '\n';
          intPoi += dstPoi;
        }
        ++bar;
      }
    }

    bar.finish();
    if (verbose) std::cout << "done (" << timer.elapsed() << "ms)." << std::endl;
  }

 private:
  const GraphT& graph; // The network we work on.
  double totPop;       // The total number of inhabitants living in the network.
  double totPoi;       // The total number of opportunities in the network.
  const int seed;      // The seed with which the random number generator will be started.
  const bool verbose;  // Should we display informative messages?
};
