#pragma once

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
#include "Tools/OpenMP.h"
#include "Tools/Timer.h"

// A travel demand calculator based on evaluating radiation model's formula for each pair of
// vertices. For each vertex s in the graph, we run a Dijkstra search without early stopping.
// Whenever we settle a vertex u, we know the intervening population between s and u, and we
// evaluate the formula for s and u.
template <typename GraphT>
class FormulaDemandCalculator {
 public:
  // Constructs a travel demand calculator for the specified network.
  explicit FormulaDemandCalculator(const GraphT& graph, const bool verbose = false) noexcept
      : graph(graph), totPop(0), verbose(verbose) {
    FORALL_VERTICES(graph, v)
      totPop += graph.population(v);
    assert(totPop > 0);
  }

  // Generates OD pairs and writes them to the specified file.
  void calculateDemand(int numODPairs, double lambda, const std::string& fileName) const {
    Timer timer;
    if (verbose) std::cout << "Calculating demand: ";
    ProgressBar bar(graph.numVertices(), verbose);

    auto gamma = numODPairs / totPop;
    assert(gamma >= 0); assert(gamma <= 1);
    assert(lambda >= 0); assert(lambda <= 1);

    using LabelSet = BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>;
    using Dijkstra = StandardDijkstra<GraphT, TravelTimeAttribute, LabelSet>;

    #pragma omp parallel
    {
      Dijkstra dijkstra(graph);
      std::minstd_rand rand(omp_get_thread_num() + 1);

      std::ofstream out(fileName + ".part" + std::to_string(omp_get_thread_num()));
      assert(out.good());

      #pragma omp for schedule(static, 1) nowait
      for (auto src = 0; src < graph.numVertices(); ++src) {
        auto srcPop = graph.population(src); // The source population.
        auto intPop = 0;                     // The intervening population.

        auto normConst = 1 / (1 -
            (1 - std::pow(lambda, totPop)) / (1 - std::pow(lambda, srcPop)) * (srcPop / totPop));
        auto m = (1 - std::pow(lambda, srcPop)) / srcPop;

        dijkstra.init({{src}});
        dijkstra.settleNextVertex();
        while (!dijkstra.queue.empty()) {
          auto dst = dijkstra.settleNextVertex();
          auto dstPop = graph.population(dst);
          if (dstPop == 0)
            continue;

          // Evaluate radiation model's formula for src and dst.
          auto ms = (1 - std::pow(lambda, srcPop + intPop)) / (srcPop + intPop);
          auto mms = (1 - std::pow(lambda, srcPop + dstPop + intPop)) / (srcPop + dstPop + intPop);
          auto prob = (ms - mms) / m;
          auto numTravelers = std::binomial_distribution<>(srcPop, gamma * normConst * prob)(rand);

          for (auto i = 0; i < numTravelers; ++i)
            out << src << ',' << dst << '\n';
          intPop += dstPop;
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
  const bool verbose;  // Should we display informative messages?
};
