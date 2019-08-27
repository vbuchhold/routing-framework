#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <csv.h>

#include "DataStructures/Graph/Attributes/CoordinateAttribute.h"
#include "DataStructures/Graph/Attributes/SequentialVertexIdAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "RawData/Visum/ZonePolygons.h"
#include "Tools/CommandLine/CommandLineParser.h"
#include "Tools/CommandLine/ProgressBar.h"
#include "Tools/Constants.h"
#include "Tools/DateHelpers.h"
#include "Tools/LexicalCast.h"
#include "Tools/StringHelpers.h"

inline void printUsage() {
  std::cout <<
      "Usage: ParseMobiToppFile -g <file> -v <file> -m <file> -o <file>\n"
      "Extracts all OD pairs from a mobiTopp file that utilize a given transportation\n"
      "mode and fall within a given analysis period.\n"
      "  -p <factor>       multiply coordinates in files by <factor> (defaults to 1)\n"
      "  -s <seed>         start random number generator with <seed> (defaults to 0)\n"
      "  -tm <mode>        extract OD pairs that utilize mode (defaults to MIV)\n"
      "  -ed <day>         earliest departure day such as Tue\n"
      "  -et <hh:mm:ss>    earliest departure time such as 07:00:00\n"
      "  -ld <day>         latest departure day such as Tue\n"
      "  -lt <hh:mm:ss>    latest departure time such as 09:00:00\n"
      "  -col <name>       apply filter to column <name> in zone file\n"
      "  -val <values>     space-separated list of values to be filtered for\n"
      "  -g <file>         network of interest in binary format\n"
      "  -v <file>         directory containing the Visum network files\n"
      "  -m <file>         mobiTopp file containing the OD pairs\n"
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
    const auto requiredMode = clp.getValue<std::string>("tm", "MIV");
    const auto earliestDay = clp.getValue<std::string>("ed");
    const auto latestDay = clp.getValue<std::string>("ld");
    const auto column = clp.getValue<std::string>("col", "CODE");
    const auto vals = clp.getValues<std::string>("val");
    const auto graphFileName = clp.getValue<std::string>("g");
    const auto visumFileName = clp.getValue<std::string>("v");
    const auto mobiToppFileName = clp.getValue<std::string>("m");
    auto outputFileName = clp.getValue<std::string>("o");
    if (!endsWith(outputFileName, ".csv"))
      outputFileName += ".csv";

    auto earliestTime = clp.getValue<std::string>("et");
    auto latestTime = clp.getValue<std::string>("lt");
    const auto earliestDep = secondsSinceMonMidnight(earliestDay, earliestTime);
    const auto latestDep = secondsSinceMonMidnight(latestDay, latestTime);

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

    // Write header to the output file.
    std::ofstream outputFile(outputFileName);
    if (!outputFile.good())
      throw std::invalid_argument("file cannot be opened -- '" + outputFileName + "'");
    outputFile << "# Input graph: " << graphFileName << '\n';
    outputFile << "# Methodology: mobiTopp\n";
    outputFile << "origin,destination,origin_zone,destination_zone,departure\n";

    // Extract the OD pairs.
    std::cout << "Extracting OD pairs..." << std::flush;
    std::minstd_rand rand(seed + 1);
    int householdId, agentId, travelTime, purpose, durOfActivity;
    char *dayOfWeek = nullptr, *dep = nullptr, *mode = nullptr;
    char *arr, *srcZone, *dstZone, *length;
    io::CSVReader<12, io::trim_chars<>, io::no_quote_escape<';'>> mobiToppFile(mobiToppFileName);
    mobiToppFile.set_header(
        "household_id", "agent_id", "day_of_week", "departure_time", "arrival_time", "travel_time",
        "origin_zone", "destination_zone", "length", "mode", "purpose", "duration_of_activity");
    while (mobiToppFile.read_row(householdId, agentId, dayOfWeek, dep, arr, travelTime,
        srcZone, dstZone, length, mode, purpose, durOfActivity)) {
      // Does the current OD pair fall within the required analysis period?
      const auto depart = secondsSinceMonMidnight(dayOfWeek, dep);
      if (depart < earliestDep || latestDep <= depart)
        continue;

      // Does the current OD pair utilize the required transportation mode?
      if (std::strcmp(mode, requiredMode.c_str()))
        continue;

      // Does the current OD pair fall within the area of interest?
      assert(srcZone[0] == 'Z');
      assert(dstZone[0] == 'Z');
      const auto srcZoneId = lexicalCast<int>(++srcZone);
      const auto dstZoneId = lexicalCast<int>(++dstZone);
      const auto srcRange = vertexRanges.find(srcZoneId);
      const auto dstRange = vertexRanges.find(dstZoneId);
      if (srcRange == vertexRanges.end() || srcRange->second.first == srcRange->second.second ||
          dstRange == vertexRanges.end() || dstRange->second.first == dstRange->second.second)
        continue;

      // Write the current OD pair to the output file.
      std::uniform_int_distribution<> srcDist(srcRange->second.first, srcRange->second.second - 1);
      std::uniform_int_distribution<> dstDist(dstRange->second.first, dstRange->second.second - 1);
      const auto s = verticesByZone[srcDist(rand)];
      const auto t = verticesByZone[dstDist(rand)];
      outputFile << s << ',' << t << ',' << srcZoneId << ',' << dstZoneId << ',' << depart << '\n';
    }
  std::cout << " done.\n";
  } catch (std::exception& e) {
    std::cerr << argv[0] << ": " << e.what() << '\n';
    std::cerr << "Try '" << argv[0] <<" -help' for more information.\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
