#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include <csv.h>

#include "Algorithms/DemandCalculation/ChooserDemandCalculator.h"
#include "Algorithms/DemandCalculation/DijkstraOpportunityChooser.h"
#include "Algorithms/DemandCalculation/FormulaDemandCalculator.h"
#include "Algorithms/DemandCalculation/KDTreeOpportunityChooser.h"
#include "Algorithms/DemandCalculation/PopulationAssignment.h"
#include "DataStructures/Graph/Attributes/LatLngAttribute.h"
#include "DataStructures/Graph/Attributes/PopulationAttribute.h"
#include "DataStructures/Graph/Attributes/SequentialVertexIdAttribute.h"
#include "DataStructures/Graph/Attributes/TravelTimeAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "Tools/CommandLine/CommandLineParser.h"
#include "Tools/Constants.h"
#include "Tools/StringHelpers.h"

inline void printUsage() {
  std::cout <<
      "Usage: CalculateDemand -n <num> -f <fmt> -g <file> -p <file> -o <file>\n"
      "Calculates travel demand in a road network according to the radiation model with\n"
      "selection, which requires as input only a population grid.\n"
      "  -n <num>          number of OD pairs to be generated\n"
      "  -l <real>         radiation model's parameter lambda (defaults to 0.999988)\n"
      "  -r <num>          use Moore neighborhoods of max range <num> (defaults to 1)\n"
      "  -f <fmt>          format of the population grid file\n"
      "                      possible values: DE EU\n"
      "  -a <algo>         use algorithm <algo> to calculate travel demand\n"
      "                      possible values: formula Dij (default) kd-tree\n"
      "  -g <file>         network of interest\n"
      "  -p <file>         population grid\n"
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

    // Parse the command-line options.
    const auto numODPairs = clp.getValue<int>("n");
    const auto lambda = clp.getValue<double>("l", 0.999988);
    const auto maxRange = clp.getValue<int>("r", 1);
    const auto gridFileFormat = clp.getValue<std::string>("f");
    const auto algo = clp.getValue<std::string>("a", "Dij");
    const auto graphFileName = clp.getValue<std::string>("g");
    const auto gridFileName = clp.getValue<std::string>("p");
    auto outputFileName = clp.getValue<std::string>("o");
    if (!endsWith(outputFileName, ".csv"))
      outputFileName += ".csv";
    const auto partFileStem = "/tmp/" + outputFileName.substr(outputFileName.rfind('/') + 1);

    // Validate command-line options.
    if (numODPairs <= 0)
      throw std::invalid_argument("number of OD pairs is no larger than 0");
    if (lambda < 0)
      throw std::invalid_argument("lambda is smaller than 0");
    if (lambda >= 1)
      throw std::invalid_argument("lambda is no smaller than 1");
    if (maxRange < 0)
      throw std::invalid_argument("max range is smaller than 0");

    // Read the graph from file.
    std::cout << "Reading graph from file..." << std::flush;
    using VertexAttributes = VertexAttrs<
        LatLngAttribute, PopulationAttribute, SequentialVertexIdAttribute>;
    using EdgeAttributes = EdgeAttrs<TravelTimeAttribute>;
    using GraphT = StaticGraph<VertexAttributes, EdgeAttributes>;
    std::ifstream graphFile(graphFileName, std::ios::binary);
    if (!graphFile.good())
      throw std::invalid_argument("file not found -- '" + graphFileName + "'");
    GraphT graph(graphFile);
    graphFile.close();
    if (graph.numVertices() > 0 && graph.sequentialVertexId(0) == INVALID_VERTEX)
      FORALL_VERTICES(graph, v) {
        assert(graph.sequentialVertexId(v) == INVALID_VERTEX);
        graph.sequentialVertexId(v) = v;
      }
    std::cout << " done.\n";

    // Assign the population grid to the graph.
    if (gridFileFormat == "DE") {
      using GridFileReaderT = io::CSVReader<2, io::trim_chars<>, io::no_quote_escape<';'>>;
      GridFileReaderT gridFileReader(gridFileName);
      gridFileReader.read_header(io::ignore_extra_column, "Gitter_ID_100m", "Einwohner");
      PopulationAssignment<GraphT, GridFileReaderT> assign(graph, gridFileReader, 100, maxRange);
      assign.run(true);
    } else if (gridFileFormat == "EU") {
      using GridFileReaderT = io::CSVReader<2>;
      GridFileReaderT gridFileReader(gridFileName);
      gridFileReader.read_header(io::ignore_extra_column, "GRD_ID", "TOT_P");
      PopulationAssignment<GraphT, GridFileReaderT> assign(graph, gridFileReader, 1000, maxRange);
      assign.run(true);
    } else {
      throw std::invalid_argument("invalid grid file format -- '" + gridFileFormat + "'");
    }

    // Calculate travel demand using the specified algorithm.
    if (algo == "formula") {
      FormulaDemandCalculator<GraphT> calculator(graph, true);
      calculator.calculateDemand(numODPairs, lambda, partFileStem);
    } else if (algo == "Dij") {
      ChooserDemandCalculator<GraphT, DijkstraOpportunityChooser> calculator(graph, true);
      calculator.calculateDemand(numODPairs, lambda, partFileStem);
    } else if (algo == "kd-tree") {
      ChooserDemandCalculator<GraphT, KDTreeOpportunityChooser> calculator(graph, true);
      calculator.calculateDemand(numODPairs, lambda, partFileStem);
    } else {
      throw std::invalid_argument("invalid algorithm -- '" + algo + "'");
    }

    // Merge the part files into a single output file.
    std::cout << "Merging part files into single output file..." << std::flush;
    std::ofstream outputFile(outputFileName);
    if (!outputFile.good())
      throw std::invalid_argument("file cannot be opened -- '" + outputFileName + "'");
    outputFile << "# Input graph: " << graphFileName << '\n';
    outputFile << "# Methodology: radiation model with selection (";
    outputFile << algo << ", " << gridFileFormat << ", l=" << lambda << ", r=" << maxRange << ")\n";
    outputFile << "origin,destination\n";
    int src, dst;
    for (auto i = 0; true; ++i) {
      const auto partFileName = partFileStem + ".part" + std::to_string(i);
      std::ifstream partFile(partFileName);
      if (!partFile.good())
        break;
      io::CSVReader<2> partFileReader(partFileName, partFile);
      partFileReader.set_header("origin", "destination");
      while (partFileReader.read_row(src, dst))
        outputFile << graph.sequentialVertexId(src) << ',' << graph.sequentialVertexId(dst) << '\n';
      partFile.close();
      std::remove(partFileName.c_str());
    }
    std::cout << " done.\n";
  } catch (std::exception& e) {
    std::cerr << argv[0] << ": " << e.what() << '\n';
    std::cerr << "Try '" << argv[0] <<" -help' for more information.\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
