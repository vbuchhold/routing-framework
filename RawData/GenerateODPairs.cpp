#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

#include "DataStructures/Graph/Attributes/LengthAttribute.h"
#include "DataStructures/Graph/Attributes/SequentialVertexIdAttribute.h"
#include "DataStructures/Graph/Attributes/TravelTimeAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Utilities/OriginDestination.h"
#include "Experiments/ODPairGenerator.h"
#include "Tools/CommandLine/CommandLineParser.h"
#include "Tools/CommandLine/ProgressBar.h"
#include "Tools/Timer.h"

void printUsage() {
  std::cout <<
      "Usage: GenerateODPairs -n <num> -i <file> -o <file>\n"
      "       GenerateODPairs -n <num> -i <file> -o <file> -r <rank> [-geom] [-len]\n"
      "       GenerateODPairs -n <num> -i <file> -o <file> -d <dist> [-geom] [-len]\n"
      "Generate OD-pairs with the origin chosen uniformly at random. The destination\n"
      "is also picked uniformly at random, or chosen by distance or Dijkstra rank.\n"
      "  -len              use physical length as cost function (default: travel time)\n"
      "  -n <num>          generate <num> OD-pairs (per rank/distance)\n"
      "  -s <seed>         start the random number generator with <seed>\n"
      "  -r <rank>         (expected) Dijkstra rank\n"
      "  -d <dist>         (expected) distance between a pair's origin and destination\n"
      "  -geom             geometrically distributed ranks/distances\n"
      "  -i <file>         the input graph in binary format\n"
      "  -o <file>         place output in <file>\n"
      "  -help             display this help and exit\n";
}

int main(int argc, char* argv[]) {
  try {
    CommandLineParser clp(argc, argv);
    if (clp.isSet("help")) {
      printUsage();
      return EXIT_SUCCESS;
    }

    std::cout << "Reading input graph from file..." << std::flush;
    const auto infile = clp.getValue<std::string>("i");
    using VertexAttributes = VertexAttrs<SequentialVertexIdAttribute>;
    using EdgeAttributes = EdgeAttrs<LengthAttribute, TravelTimeAttribute>;
    using Graph = StaticGraph<VertexAttributes, EdgeAttributes>;
    std::ifstream in(infile, std::ios::binary);
    if (!in.good())
      throw std::invalid_argument("file not found -- '" + infile + "'");
    Graph graph(in);
    in.close();
    std::cout << " done." << std::endl;
    std::default_random_engine rand(clp.getValue<int>("s", 19900325));
    ODPairGenerator<Graph, TravelTimeAttribute> gen(graph, rand);

    if (graph.numVertices() > 0 && graph.sequentialVertexId(0) == INVALID_VERTEX)
      FORALL_VERTICES(graph, v)
        graph.sequentialVertexId(v) = v;

    const bool useLen = clp.isSet("len");
    if (useLen)
      // Use physical length as cost function.
      FORALL_EDGES(graph, e)
        graph.travelTime(e) = graph.length(e);

    const auto outfile = clp.getValue<std::string>("o");
    std::ofstream out(outfile + ".csv");
    if (!out.good())
      throw std::invalid_argument("file cannot be opened -- '" + outfile + ".csv'");
    out << "# Input graph: " << infile << std::endl;
    out << "# Methodology: ";

    const int numPairs = clp.getValue<int>("n");
    ProgressBar bar;
    bar.setPercentageOutputInterval(25);

    if (clp.isSet("r")) {
      // Choose the destination by Dijkstra rank.
      const auto ranks = clp.getValues<int>("r");
      const bool isGeom = clp.isSet("geom");

      std::cout << "The destinations are chosen by Dijkstra rank (";
      std::cout << (isGeom ? "geometrically distributed" : "constant") << ")." << std::endl;
      std::cout << "Cost function: " << (useLen ? "physical length" : "travel time") << std::endl;

      out << (isGeom ? "geometrically distributed" : "constant") << " Dijkstra ranks (";
      for (int i = 0; i < ranks.size(); ++i) {
        if (ranks[i] < 0)
          throw std::invalid_argument("negative exponent -- '" + std::to_string(ranks[i]) + "'");
        out << (i > 0 ? " " : "") << ranks[i];
      }
      out << ")\n";
      out << "origin,destination,dijkstra_rank\n";

      std::vector<std::pair<OriginDestination, int>> result;
      result.reserve(numPairs * ranks.size());
      Timer timer;
      for (const int r : ranks) {
        std::cout << "Generating " << numPairs << " OD-pairs (" << r << "): ";
        bar.init(numPairs);
        std::geometric_distribution<> dist(1.0 / r);
        for (int i = 0; i < numPairs; ++i) {
          const int actualRank = std::min(isGeom ? dist(rand) : r, graph.numVertices() - 1);
          result.emplace_back(gen.getRandomODPairChosenByDijkstraRank(actualRank), actualRank);
          ++bar;
        }
        std::cout << "done." << std::endl;
      }
      const int elapsed = timer.elapsed();
      std::cout << "Running time: " << elapsed << "ms" << std::endl;

      std::cout << "Writing OD-pairs to file..." << std::flush;
      for (const auto& record : result) {
        out << graph.sequentialVertexId(record.first.origin) << ',';
        out << graph.sequentialVertexId(record.first.destination) << ',';
        out << record.second << '\n';
      }
      std::cout << " done." << std::endl;
    } else if (clp.isSet("d")) {
      // Choose the destination by distance.
      const auto distances = clp.getValues<int>("d");
      const bool isGeom = clp.isSet("geom");

      std::cout << "The destinations are chosen by distance (";
      std::cout << (isGeom ? "geometrically distributed" : "constant") << ")." << std::endl;
      std::cout << "Cost function: " << (useLen ? "physical length" : "travel time") << std::endl;

      out << (isGeom ? "geometrically distributed" : "constant") << " distances (";
      for (int i = 0; i < distances.size(); ++i) {
        const int d = distances[i];
        if (d < 0)
          throw std::invalid_argument("negative distance -- '" + std::to_string(d) + "'");
        out << (i > 0 ? " " : "") << d;
      }
      out << ")\n";
      out << "origin,destination\n";

      std::vector<OriginDestination> result;
      result.reserve(numPairs * distances.size());
      Timer timer;
      for (const int d : distances) {
        std::cout << "Generating " << numPairs << " OD-pairs (" << d << "): ";
        bar.init(numPairs);
        std::geometric_distribution<> dist(1.0 / d);
        for (int i = 0; i < numPairs; ++i) {
          result.push_back(gen.getRandomODPairChosenByDistance(isGeom ? dist(rand) : d));
          ++bar;
        }
        std::cout << "done." << std::endl;
      }
      const int elapsed = timer.elapsed();
      std::cout << "Running time: " << elapsed << "ms" << std::endl;

      std::cout << "Writing OD-pairs to file..." << std::flush;
      for (const auto& record : result) {
        out << graph.sequentialVertexId(record.origin) << ',';
        out << graph.sequentialVertexId(record.destination) << '\n';
      }
      std::cout << " done." << std::endl;
    } else {
      // Choose the destination uniformly at random.
      std::cout << "The destinations are chosen uniformly at random." << std::endl;
      out << "random\n";
      out << "origin,destination\n";

      std::vector<OriginDestination> result;
      result.reserve(numPairs);
      Timer timer;
      std::cout << "Generating " << numPairs << " OD-pairs: ";
      bar.init(numPairs);
      for (int i = 0; i < numPairs; ++i) {
        result.push_back(gen.getRandomODPair());
        ++bar;
      }
      std::cout << "done." << std::endl;
      const int elapsed = timer.elapsed();
      std::cout << "Running time: " << elapsed << "ms" << std::endl;

      std::cout << "Writing OD-pairs to file..." << std::flush;
      for (const auto& record : result) {
        out << graph.sequentialVertexId(record.origin) << ',';
        out << graph.sequentialVertexId(record.destination) << '\n';
      }
      std::cout << " done." << std::endl;
    }
  } catch (std::exception& e) {
    std::cerr << argv[0] << ": " << e.what() << std::endl;
    std::cerr << "Try '" << argv[0] <<" -help' for more information." << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
