#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "DataStructures/Graph/Attributes/CoordinateAttribute.h"
#include "DataStructures/Graph/Attributes/SequentialVertexIdAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "RawData/Visum/ZonePolygons.h"
#include "Tools/CommandLine/CommandLineParser.h"
#include "Tools/CommandLine/ProgressBar.h"
#include "Tools/Constants.h"
#include "Tools/StringHelpers.h"
#include "Tools/Workarounds.h"

inline void printUsage() {
  std::cout <<
      "Usage: ParseMtxFile [-n <num>] -g <file> -v <file> -m <file> -o <file>\n"
      "Reads a MTX file containing an OD matrix with fractional demands and samples a\n"
      "set of OD pairs that have the same distribution.\n"
      "  -n <num>          number of OD pairs to be sampled (default: take from MTX)\n"
      "  -p <factor>       multiply coordinates in files by <factor> (defaults to 1)\n"
      "  -s <seed>         start random number generator with <seed> (defaults to 0)\n"
      "  -col <name>       apply filter to column <name> in zone file\n"
      "  -val <values>     space-separated list of values to be filtered for\n"
      "  -g <file>         network of interest in binary format\n"
      "  -v <file>         directory containing the Visum network files\n"
      "  -m <file>         MTX file containing the OD matrix\n"
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
    const auto precision = clp.getValue<int>("p", 1);
    const auto seed = clp.getValue<int>("s", 0);
    const auto column = clp.getValue<std::string>("col", "CODE");
    const auto vals = clp.getValues<std::string>("val");
    const auto graphFileName = clp.getValue<std::string>("g");
    const auto visumFileName = clp.getValue<std::string>("v");
    const auto mtxFileName = clp.getValue<std::string>("m");
    auto outputFileName = clp.getValue<std::string>("o");
    if (!endsWith(outputFileName, ".csv"))
      outputFileName += ".csv";

    // Read graph and zone polygons from file.
    std::cout << "Reading graph..." << std::flush;
    std::ifstream graphFile(graphFileName, std::ios::binary);
    if (!graphFile.good())
      throw std::invalid_argument("file not found -- '" + graphFileName + "'");
    StaticGraph<VertexAttrs<CoordinateAttribute, SequentialVertexIdAttribute>> graph(graphFile);
    graphFile.close();
    if (graph.numVertices() > 0 && graph.sequentialVertexId(0) == INVALID_VERTEX)
      FORALL_VERTICES(graph, v) {
        assert(graph.sequentialVertexId(v) == INVALID_VERTEX);
        graph.sequentialVertexId(v) = v;
      }
    std::cout << " done.\n";
    const auto zonePolygons = visum::readZonePolygonsFrom(visumFileName, precision, column, vals);

    // Assign vertices to zones.
    std::cout << "Assigning vertices to zones: ";
    ProgressBar bar(zonePolygons.size(), true);
    std::unordered_map<int, std::pair<int, int>> vertexRanges;
    std::vector<int> verticesByZone;
    for (const auto& zone : zonePolygons) {
      auto& range = vertexRanges[zone.first];
      range.first = verticesByZone.size();
      const auto box = zone.second.boundingBox();
      FORALL_VERTICES(graph, v)
        if (box.contains(graph.coordinate(v)) && zone.second.contains(graph.coordinate(v)))
          verticesByZone.push_back(graph.sequentialVertexId(v));
      range.second = verticesByZone.size();
      ++bar;
    }
    std::cout << "done.\n";

    std::vector<int> srcZones;
    std::vector<int64_t> totalOutflows;
    std::vector<int> sampledOutflows;
    std::vector<int> firstOutflow;
    std::vector<int> dstZones;
    std::vector<int> outflows;

    // Read the MTX file.
    std::cout << "Reading MTX file..." << std::flush;
    std::ifstream mtxFile(mtxFileName);
    if (!mtxFile.good())
      throw std::invalid_argument("file not found -- '" + mtxFileName + "'");
    for (auto i = 0; i < 8; ++i)
      mtxFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    int src, dst;
    double outflow;
    auto prevSrc = -1, prevDst = -1;
    unused(prevSrc, prevDst);
    while (mtxFile >> src >> dst >> outflow) {
      assert(src >= 0);
      assert(dst >= 0);
      assert(src >= prevSrc);
      assert(src != prevSrc || dst > prevDst);
      const auto srcRange = vertexRanges.find(src);
      const auto dstRange = vertexRanges.find(dst);
      if (srcRange != vertexRanges.end() && srcRange->second.first != srcRange->second.second &&
          dstRange != vertexRanges.end() && dstRange->second.first != dstRange->second.second &&
          outflow > 0) {
        if (srcZones.empty() || src != srcZones.back()) {
          srcZones.push_back(src);
          firstOutflow.push_back(outflows.size());
        }
        dstZones.push_back(dst);
        outflows.push_back(std::round(outflow * 1000));
      }
      prevSrc = src;
      prevDst = dst;
    }
    firstOutflow.push_back(outflows.size());
    std::cout << " done.\n";

    // Form the prefix sums of the outflows.
    std::cout << "Forming prefix sums..." << std::flush;
    totalOutflows.assign(srcZones.size(), 0);
    sampledOutflows.assign(srcZones.size(), 0);
    for (auto i = 0; i < srcZones.size(); ++i) {
      const auto first = outflows.begin() + firstOutflow[i];
      const auto last = outflows.begin() + firstOutflow[i + 1];
      std::partial_sum(first, last, first);
      totalOutflows[i] = *(last - 1);
    }
    std::partial_sum(totalOutflows.begin(), totalOutflows.end(), totalOutflows.begin());
    std::cout << " done.\n";

    // Sample the source zones.
    const auto numODPairs = clp.getValue<int>("n", std::round(totalOutflows.back() / 1000.0));
    std::cout << "Sampling source zones: ";
    bar.init(numODPairs);
    std::minstd_rand rand(seed + 1);
    std::uniform_int_distribution<int64_t> rankDist(0, totalOutflows.back() - 1);
    for (auto i = 0; i < numODPairs; ++i) {
      const auto rank = rankDist(rand);
      auto j = 0;
      while (totalOutflows[j] <= rank) ++j;
      assert(j < totalOutflows.size());
      ++sampledOutflows[j];
      ++bar;
    }
    std::cout << "done.\n";

    // Sample the OD pairs.
    std::cout << "Sampling OD pairs: ";
    bar.init(numODPairs);
    std::ofstream outputFile(outputFileName);
    if (!outputFile.good())
      throw std::invalid_argument("file cannot be opened -- '" + outputFileName + "'");
    outputFile << "# Input graph: " << graphFileName << '\n';
    outputFile << "# Methodology: Visum\n";
    outputFile << "origin,destination,origin_zone,destination_zone\n";
    for (auto i = 0; i < srcZones.size(); ++i) {
      const auto srcZone = srcZones[i];
      const auto& srcRange = vertexRanges[srcZone];
      assert(srcRange.first < srcRange.second);
      std::uniform_int_distribution<> srcDist(srcRange.first, srcRange.second - 1);
      std::uniform_int_distribution<> rankDist(0, outflows[firstOutflow[i + 1] - 1] - 1);

      for (auto k = 0; k < sampledOutflows[i]; ++k) {
        const auto rank = rankDist(rand);
        auto j = firstOutflow[i];
        while (outflows[j] <= rank) ++j;
        assert(j < firstOutflow[i + 1]);

        const auto dstZone = dstZones[j];
        const auto& dstRange = vertexRanges[dstZone];
        assert(dstRange.first < dstRange.second);
        std::uniform_int_distribution<> dstDist(dstRange.first, dstRange.second - 1);

        const auto srcVertex = verticesByZone[srcDist(rand)];
        const auto dstVertex = verticesByZone[dstDist(rand)];
        outputFile << srcVertex << ',' << dstVertex << ',' << srcZone << ',' << dstZone << '\n';
        ++bar;
      }
    }
    std::cout << "done.\n";
  } catch (std::exception& e) {
    std::cerr << argv[0] << ": " << e.what() << '\n';
    std::cerr << "Try '" << argv[0] <<" -help' for more information.\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
