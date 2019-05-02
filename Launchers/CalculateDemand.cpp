#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <csv.h>

#include "Algorithms/DemandCalculation/ChooserDemandCalculator.h"
#include "Algorithms/DemandCalculation/DijkstraOpportunityChooser.h"
#include "Algorithms/DemandCalculation/FormulaDemandCalculator.h"
#include "Algorithms/DemandCalculation/KDTreeOpportunityChooser.h"
#include "Algorithms/DemandCalculation/PopulationAssignment.h"
#include "DataStructures/Containers/BitVector.h"
#include "DataStructures/Geometry/Area.h"
#include "DataStructures/Geometry/KDTree.h"
#include "DataStructures/Geometry/Point.h"
#include "DataStructures/Graph/Attributes/LatLngAttribute.h"
#include "DataStructures/Graph/Attributes/LengthAttribute.h"
#include "DataStructures/Graph/Attributes/NumOpportunitiesAttribute.h"
#include "DataStructures/Graph/Attributes/PopulationAttribute.h"
#include "DataStructures/Graph/Attributes/SequentialVertexIdAttribute.h"
#include "DataStructures/Graph/Attributes/TravelTimeAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "Tools/CommandLine/CommandLineParser.h"
#include "Tools/Constants.h"
#include "Tools/StringHelpers.h"

inline void printUsage() {
  std::cout <<
      "Usage: CalculateDemand -n <num> -f <fmt> -g <file> -pop <file> -o <file>\n"
      "Calculates travel demand in a road network according to the radiation model with\n"
      "selection, which requires as input only a population grid.\n"
      "  -len              use physical lengths as metric (default: travel time)\n"
      "  -n <num>          number of OD pairs to be generated\n"
      "  -l <real>         radiation model's parameter lambda (defaults to 0.999988)\n"
      "  -p <prob>         swap src and dst of each OD pair with <prob> (defaults to 0)\n"
      "  -r <num>          use Moore neighborhoods of max range <num> (defaults to 1)\n"
      "  -d <meters>       do not assign POIs to vertices farther than <meters>\n"
      "  -s <seed>         start demand calculation with <seed> (defaults to 0)\n"
      "  -f <fmt>          format of the population grid file\n"
      "                      possible values: DE EU\n"
      "  -m <algo>         use method <algo> to calculate travel demand\n"
      "                      possible values: formula Dij (default) kd-tree\n"
      "  -g <file>         input graph in binary format\n"
      "  -a <file>         restrict origins and destinations to polygonal study area\n"
      "  -pop <file>       population grid\n"
      "  -poi <file>       use density of POIs as proxy for trip attraction rates\n"
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
    const auto useLengths = clp.isSet("len");
    const auto numODPairs = clp.getValue<int>("n");
    const auto lambda = clp.getValue<double>("l", 0.999988);
    const auto swapProb = clp.getValue<double>("p", 0.0);
    const auto maxRange = clp.getValue<int>("r", 1);
    const auto maxDistance = clp.getValue<int>("d", 200);
    const auto seed = clp.getValue<int>("s", 0);
    const auto popFileFormat = clp.getValue<std::string>("f");
    const auto algo = clp.getValue<std::string>("m", "Dij");
    const auto graphFileName = clp.getValue<std::string>("g");
    const auto areaFileName = clp.getValue<std::string>("a");
    const auto popFileName = clp.getValue<std::string>("pop");
    const auto poiFileName = clp.getValue<std::string>("poi");
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
    if (swapProb > 1)
      throw std::invalid_argument("swap probability is larger than 1");
    if (maxRange < 0)
      throw std::invalid_argument("max range is smaller than 0");
    if (maxDistance < 0)
      throw std::invalid_argument("max distance is smaller than 0");
    if (seed < 0)
      throw std::invalid_argument("seed is smaller than 0");

    // Read the graph from file.
    std::cout << "Reading graph from file..." << std::flush;
    using VertexAttributes = VertexAttrs<
        LatLngAttribute, NumOpportunitiesAttribute,
        PopulationAttribute, SequentialVertexIdAttribute>;
    using EdgeAttributes = EdgeAttrs<LengthAttribute, TravelTimeAttribute>;
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
    if (useLengths)
      FORALL_EDGES(graph, e)
        graph.travelTime(e) = graph.length(e);
    std::cout << " done.\n";

    // Read the study area from OSM POLY file.
    BitVector isVertexInStudyArea(graph.numVertices(), true);
    if (!areaFileName.empty()) {
      std::cout << "Reading study area from OSM POLY file..." << std::flush;
      Area studyArea;
      studyArea.importFromOsmPolyFile(areaFileName);
      const auto box = studyArea.boundingBox();
      FORALL_VERTICES(graph, u) {
        const Point p(graph.latLng(u).longitude(), graph.latLng(u).latitude());
        isVertexInStudyArea[u] = box.contains(p) && studyArea.contains(p);
      }
      std::cout << " done.\n";
    }

    // Assign the population grid to the graph.
    if (popFileFormat == "DE") {
      using PopFileReaderT = io::CSVReader<2, io::trim_chars<>, io::no_quote_escape<';'>>;
      using PopAssignmentT = PopulationAssignment<GraphT, PopFileReaderT>;
      PopFileReaderT popFileReader(popFileName);
      popFileReader.read_header(io::ignore_extra_column, "Gitter_ID_100m", "Einwohner");
      PopAssignmentT assign(graph, isVertexInStudyArea, popFileReader, 100, maxRange);
      assign.run(true);
    } else if (popFileFormat == "EU") {
      using PopFileReaderT = io::CSVReader<2>;
      using PopAssignmentT = PopulationAssignment<GraphT, PopFileReaderT>;
      PopFileReaderT popFileReader(popFileName);
      popFileReader.read_header(io::ignore_extra_column, "GRD_ID", "TOT_P");
      PopAssignmentT assign(graph, isVertexInStudyArea, popFileReader, 1000, maxRange);
      assign.run(true);
    } else {
      throw std::invalid_argument("invalid population file format -- '" + popFileFormat + "'");
    }

    // Assign POIs to vertices (or use the density of population as a proxy for trip attraction).
    if (!poiFileName.empty()) {
      Timer timer;
      std::cout << "Assigning POIs to vertices..." << std::flush;
      int assignedPoi = 0;
      int unassignedPoi = 0;

      std::vector<int> verticesInStudyArea(isVertexInStudyArea.cardinality());
      for (auto i = 0, v = isVertexInStudyArea.firstSetBit(); i < verticesInStudyArea.size(); ++i) {
        verticesInStudyArea[i] = v;
        v = isVertexInStudyArea.nextSetBit(v);
      }

      std::vector<Point> points(verticesInStudyArea.size());
      for (auto i = 0; i < verticesInStudyArea.size(); ++i) {
        const auto& latLng = graph.latLng(verticesInStudyArea[i]);
        points[i] = {latLng.longitude(), latLng.latitude()};
      }
      KDTree tree(points);

      double lng, lat;
      using TrimPolicy = io::trim_chars<>;
      using QuotePolicy = io::no_quote_escape<','>;
      using OverflowPolicy = io::throw_on_overflow;
      using CommentPolicy = io::single_line_comment<'#'>;
      using PoiReader = io::CSVReader<2, TrimPolicy, QuotePolicy, OverflowPolicy, CommentPolicy>;
      PoiReader poiReader(poiFileName);
      poiReader.read_header(io::ignore_extra_column, "longitude", "latitude");
      while (poiReader.read_row(lng, lat)) {
        const LatLng query(lat, lng);
        const auto p = tree.findClosestPoint(Point(query.longitude(), query.latitude()));
        if (graph.latLng(verticesInStudyArea[p]).getGreatCircleDistanceTo(query) <= maxDistance) {
          ++graph.numOpportunities(verticesInStudyArea[p]);
          ++assignedPoi;
        } else {
          ++unassignedPoi;
        }
      }

      std::cout << " done (" << timer.elapsed() << "ms).\n";
      std::cout << "  POIs that were assigned: " << assignedPoi << "\n";
      std::cout << "  POIs that could not be assigned: " << unassignedPoi << "\n";
    } else {
      FORALL_VERTICES(graph, v)
        graph.numOpportunities(v) = graph.population(v);
    }

    // Calculate travel demand using the specified method.
    if (algo == "formula") {
      FormulaDemandCalculator<GraphT> calculator(graph, seed, true);
      calculator.calculateDemand(numODPairs, lambda, swapProb, partFileStem);
    } else if (algo == "Dij") {
      ChooserDemandCalculator<GraphT, DijkstraOpportunityChooser> calculator(graph, seed, true);
      calculator.calculateDemand(numODPairs, lambda, swapProb, partFileStem);
    } else if (algo == "kd-tree") {
      ChooserDemandCalculator<GraphT, KDTreeOpportunityChooser> calculator(graph, seed, true);
      calculator.calculateDemand(numODPairs, lambda, swapProb, partFileStem);
    } else {
      throw std::invalid_argument("invalid method -- '" + algo + "'");
    }

    // Merge the part files into a single output file.
    std::cout << "Merging part files into single output file..." << std::flush;
    std::ofstream outputFile(outputFileName);
    if (!outputFile.good())
      throw std::invalid_argument("file cannot be opened -- '" + outputFileName + "'");
    outputFile << "# Input graph: " << graphFileName << '\n';
    outputFile << "# Methodology: radiation model with selection (";
    outputFile << algo << ", " << popFileFormat << ", l=" << lambda << ", r=" << maxRange << ")\n";
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
