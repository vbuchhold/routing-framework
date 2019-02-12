#pragma once

#include <cassert>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "DataStructures/Graph/Graph.h"
#include "Tools/CommandLine/ProgressBar.h"
#include "Tools/OpenMP.h"
#include "Tools/Timer.h"

// A travel demand calculator based on repeatedly choosing random sources and corresponding targets.
// Each source is chosen with probability proportional to its population. The target is chosen to be
// the closest opportunity with sufficiently high fitness.
template <typename GraphT, template <typename> class OpportunityChooserT>
class ChooserDemandCalculator {
 public:
  // Constructs a travel demand calculator for the specified network.
  explicit ChooserDemandCalculator(const GraphT& graph, const bool verbose = false) noexcept
      : graph(graph), totalPop(0), verbose(verbose) {
    FORALL_VERTICES(graph, v)
      totalPop += graph.population(v);
    assert(totalPop > 0);
  }

  // Generates OD pairs and writes them to the specified file.
  void calculateDemand(int numODPairs, double lambda, const std::string& fileName) const {
    Timer timer;
    if (verbose) std::cout << "Calculating demand: ";
    ProgressBar bar(numODPairs, verbose);

    // We will choose the source uniformly at random from the following box.
    // The number of objects representing a particular source is proportional to its population.
    std::vector<int> boxContainingSources(totalPop);
    for (auto v = 0, i = 0; v < graph.numVertices(); ++v)
      for (auto j = 0; j < graph.population(v); ++j)
        boxContainingSources[i++] = v;

    #pragma omp parallel
    {
      OpportunityChooserT<GraphT> chooser(graph);
      std::minstd_rand rand(omp_get_thread_num() + 1);
      std::uniform_int_distribution<> rankDist(1, totalPop);
      std::uniform_int_distribution<> sourceDist(0, totalPop - 1);

      std::ofstream out(fileName + ".part" + std::to_string(omp_get_thread_num()));
      assert(out.good());

      #pragma omp for schedule(static, 1) nowait
      for (auto i = 0; i < numODPairs; ++i) {
        auto numFitOpportunities = 0;
        while (numFitOpportunities == 0) {
          const auto rank = rankDist(rand);
          numFitOpportunities = std::binomial_distribution<>(rank, 1 - lambda)(rand);
        }
        const auto src = boxContainingSources[sourceDist(rand)];
        const auto dst = chooser.findClosestOpportunityWithHighFitness(src, numFitOpportunities);
        out << src << ',' << dst << '\n';
        ++bar;
      }
    }

    bar.finish();
    if (verbose) std::cout << "done (" << timer.elapsed() << "ms)." << std::endl;
  }

 private:
  const GraphT& graph; // The network we work on.
  int totalPop;        // The total number of inhabitants living in the network.
  const bool verbose;  // Should we display informative messages?
};
