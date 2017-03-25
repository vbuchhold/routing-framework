#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

#include "DataStructures/Graph/Attributes/TravelCostAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Utilities/OriginDestination.h"
#include "Experiments/ODPairGenerator.h"
#include "Tools/CommandLine/CommandLineParser.h"
#include "Tools/CommandLine/ProgressBar.h"

void printUsage() {
  std::cout <<
      "Usage: GenerateODPairs -n <num> [-s <seed>] -i <file> -o <file>\n"
      "       GenerateODPairs -n <num> [-s <seed>] -i <file> -o <file> -r <ranks>\n"
      "       GenerateODPairs -n <num> [-s <seed>] -i <file> -o <file> -d <dist> [-geo]\n"
      "This program generates OD-pairs, with the origin chosen uniformly at random.\n"
      "The destination is also picked uniformly at random, or chosen by distance or\n"
      "Dijkstra rank. Dijkstra ranks are specified in terms of powers of two.\n"
      "  -n <num>          the number of OD-pairs to be generated (per Dijkstra rank)\n"
      "  -s <seed>         the seed for the random number generator\n"
      "  -r <ranks>        a space-separated list of Dijkstra ranks\n"
      "  -d <dist>         (expected) distance between a pair's origin and destination\n"
      "  -geo              geometrically distributed distances with expected value -d\n"
      "  -i <file>         the input graph\n"
      "  -o <file>         the output file\n"
      "  -help             display this help and exit\n";
}

// Prints the specified error message onto standard error.
void printErrorMessage(const std::string& invokedName, const std::string& msg) {
  std::cerr << invokedName << ": " << msg << std::endl;
  std::cerr << "Try '" << invokedName <<" -help' for more information." << std::endl;
}

// Writes the specified OD-pair to the output file.
void output(const OriginDestination& pair, std::ofstream& out) {
  out << "q " << pair.origin << " " << pair.destination << std::endl;
}

int main(int argc, char* argv[]) {
  CommandLineParser clp;
  try {
    clp.parse(argc, argv);
  } catch (std::invalid_argument& e) {
    printErrorMessage(argv[0], e.what());
    return EXIT_FAILURE;
  }

  if (clp.isSet("help")) {
    printUsage();
    return EXIT_SUCCESS;
  }

  std::default_random_engine rand(clp.getValue<int>("s", 19900325));
  using Graph = StaticGraph<VertexAttrs<>, EdgeAttrs<TravelCostAttribute>>;
  std::ifstream in(clp.getValue<std::string>("i"), std::ios::binary);
  std::cout << "Reading the input graph..." << std::flush;
  Graph graph(in);
  std::cout << " done." << std::endl;
  in.close();
  ODPairGenerator<Graph, TravelCostAttribute> gen(graph, rand);

  const int numPairs = clp.getValue<int>("n");
  std::ofstream out(clp.getValue<std::string>("o") + ".od");
  out << "p od " << numPairs << std::endl;

  ProgressBar bar;
  bar.setPercentageOutputInterval(25);

  if (clp.isSet("r")) {
    // Choose the destination by Dijkstra rank.
    out << "r";
    for (const auto rank : clp.getValues<int>("r"))
      out << " " << rank;
    out << std::endl;

    std::cout << "The destinations are chosen by Dijkstra rank." << std::endl;
    for (const auto rank : clp.getValues<int>("r")) {
      std::cout << "Generating " << numPairs << " OD-pairs (2^" << rank << "): " << std::flush;
      bar.init(numPairs);
      for (int i = 0; i < numPairs; ++i) {
        output(gen.getRandomODPairChosenByDijkstraRank(std::pow(2, rank)), out);
        ++bar;
      }
      std::cout << "done." << std::endl;
    }
  } else if (clp.isSet("d")) {
    // Choose the destination by distance.
    const int distance = clp.getValue<int>("d");
    const bool isGeo = clp.isSet("geo");
    std::geometric_distribution<> dist(1.0 / distance);

    out << "d " << distance << " ";
    out << (isGeo ? "geometrically-distributed" : "exact") << std::endl;

    if (isGeo) {
      std::cout << "The origin-destination distance is geometrically distributed." << std::endl;
      std::cout << "Expected distance: " << distance << std::endl;
    } else {
      std::cout << "The origin-destination distance is " << distance << "." << std::endl;
    }

    std::cout << "Generating " << numPairs << " OD-pairs: " << std::flush;
    bar.init(numPairs);
    for (int i = 0; i < numPairs; ++i) {
      output(gen.getRandomODPairChosenByDistance(isGeo ? dist(rand) : distance), out);
      ++bar;
    }
    std::cout << "done." << std::endl;
  } else {
    // Choose the destination uniformly at random.
    std::cout << "Generating " << numPairs << " OD-pairs: " << std::flush;
    bar.init(numPairs);
    for (int i = 0; i < numPairs; ++i) {
      output(gen.getRandomODPair(), out);
      ++bar;
    }
    std::cout << "done." << std::endl;
  }

  out.close();
  return EXIT_SUCCESS;
}
