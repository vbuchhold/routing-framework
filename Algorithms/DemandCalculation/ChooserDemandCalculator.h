#pragma once

#include <cassert>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

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
  explicit ChooserDemandCalculator(const GraphT& graph, const int seed, const bool verbose) noexcept
      : graph(graph), numOpportunities(0), seed(seed), verbose(verbose) {
    FORALL_VERTICES(graph, v)
    numOpportunities += graph.numOpportunities(v);
    assert(numOpportunities > 0);
    assert(seed >= 0);
  }

  // Generates OD pairs and writes them to the specified file.
  void calculateDemand(
      int numODPairs, double lambda, double swapProb, const std::string& fileName) const {
    Timer timer;
    if (verbose) std::cout << "Calculating demand: ";
    ProgressBar bar(numODPairs, verbose);

    const auto firstWeight = &graph.population(0);
    const auto lastWeight = &graph.population(graph.numVertices() - 1) + 1;

    #pragma omp parallel
    {
      OpportunityChooserT<GraphT> chooser(graph, seed);
      std::minstd_rand rand(seed + omp_get_thread_num() + 1);
      std::discrete_distribution<> sourceDist(firstWeight, lastWeight);
      std::uniform_int_distribution<> rankDist(1, numOpportunities);
      std::bernoulli_distribution swapDist(swapProb);

      std::ofstream out(fileName + ".part" + std::to_string(omp_get_thread_num()));
      assert(out.good());

      #pragma omp for schedule(static, 1) nowait
      for (auto i = 0; i < numODPairs; ++i) {
        const auto src = sourceDist(rand);
        auto dst = src;
        while (src == dst) {
          auto numFitOpportunities = 0;
          while (numFitOpportunities == 0) {
            const auto rank = rankDist(rand);
            numFitOpportunities = std::binomial_distribution<>(rank, 1 - lambda)(rand);
          }
          dst = chooser.findClosestOpportunityWithHighFitness(src, numFitOpportunities);
        }
        if (swapDist(rand))
          out << dst << ',' << src << '\n';
        else
          out << src << ',' << dst << '\n';
        ++bar;
      }
    }

    bar.finish();
    if (verbose) std::cout << "done (" << timer.elapsed() << "ms)." << std::endl;
  }

 private:
  const GraphT& graph;  // The network we work on.
  int numOpportunities; // The total number of opportunities in the network.
  const int seed;       // The seed with which the random number generators will be started.
  const bool verbose;   // Should we display informative messages?
};
