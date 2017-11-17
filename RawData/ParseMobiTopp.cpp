#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <csv.h>

#include "DataStructures/Geometry/Area.h"
#include "DataStructures/Geometry/Rectangle.h"
#include "DataStructures/Graph/Attributes/CoordinateAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "RawData/Visum/ZonePolygons.h"
#include "Tools/CommandLine/CommandLineParser.h"
#include "Tools/CommandLine/ProgressBar.h"
#include "Tools/BinaryIO.h"
#include "Tools/DateHelpers.h"
#include "Tools/LexicalCast.h"

void printUsage() {
  std::cout <<
      "Usage: ParseMobiTopp -g <file> -v <file> -o <file>\n"
      "       ParseMobiTopp [-s <seed>] -g <file> -vs <file> -od <file> -o <file>\n"
      "This program extracts all OD-pairs between traffic zones from a given mobiTopp\n"
      "file that lie in a selected analysis period and translates the zone IDs into\n"
      "vertex IDs in the corresponding graph.\n"
      "  -s <seed>         start the random number generator with <seed>\n"
      "  -m <mode>         extract OD-pairs that use <mode> (defaults to 'MIV')\n"
      "  -ed <day>         earliest day of departure such as 'Di'\n"
      "  -et <HH:mm:ss>    earliest departure time such as 07:00:00\n"
      "  -ld <day>         latest day of departure such as 'Di'\n"
      "  -lt <HH:mm:ss>    latest departure time such as 09:00:00\n"
      "  -g <file>         corresponding graph in binary format\n"
      "  -v <file>         directory containing Visum files managing the zone polygons\n"
      "  -vs <file>        binary file containing the vertex set for each traffic zone\n"
      "  -od <file>        mobiTopp file containing OD-pairs between traffic zones\n"
      "  -o <file>         place the output into <file>\n"
      "  -help             display this help and exit\n";
}

// Computes for each traffic zone which vertices it contains.
inline void computeVertexSets(const CommandLineParser& clp) {
  const std::string graphFile = clp.getValue<std::string>("g");
  const std::string visumFile = clp.getValue<std::string>("v");
  const std::string outfile = clp.getValue<std::string>("o");

  std::cout << "Reading the graph..." << std::flush;
  std::ifstream in(graphFile, std::ios::binary);
  if (!in.good())
    throw std::invalid_argument("file not found -- '" + graphFile + "'");
  StaticGraph<VertexAttrs<CoordinateAttribute>> graph(in);
  in.close();
  std::cout << " done." << std::endl;

  std::cout << "Reading the zone polygons..." << std::flush;
  const std::map<int, Area> zoneSurfaces = visum::readZonePolygonsFrom(visumFile);
  std::cout << " done." << std::endl;

  std::cout << "Computing the vertex sets for all traffic zones: " << std::flush;
  const int numZones = zoneSurfaces.size();
  ProgressBar bar(numZones);
  bar.setPercentageOutputInterval(50);
  bar.setDotOutputInterval(10);

  std::vector<std::vector<int>> zoneVertexSets;
  for (const auto& surface : zoneSurfaces) {
    const Rectangle box = surface.second.boundingBox();
    zoneVertexSets.emplace_back();
    FORALL_VERTICES(graph, v)
      if (box.contains(graph.coordinate(v)) && surface.second.contains(graph.coordinate(v)))
        zoneVertexSets.back().push_back(v);
    ++bar;
  }
  std::cout << " done." << std::endl;

  std::cout << "Writing the vertex sets..." << std::flush;
  std::ofstream out(outfile + ".vs.bin", std::ios::binary);
  if (!out.good())
    throw std::invalid_argument("file cannot be opened -- '" + outfile + ".vs.bin'");
  bio::write(out, numZones);
  for (const auto& surface : zoneSurfaces)
    bio::write(out, surface.first);
  for (const auto& vertexSet : zoneVertexSets)
    bio::write(out, vertexSet);
  std::cout << " done." << std::endl;
}

// Extracts selected OD-pairs from a mobiTopp file.
inline void extractODPairs(const CommandLineParser& clp) {
  const int seed = clp.getValue<int>("s", 19900325);
  const std::string selectedMode = clp.getValue<std::string>("m", "MIV");
  const std::string earliestDay = clp.getValue<std::string>("ed");
  const std::string latestDay = clp.getValue<std::string>("ld");
  const std::string graphFile = clp.getValue<std::string>("g");
  const std::string vertexSetFile = clp.getValue<std::string>("vs");
  const std::string odFile = clp.getValue<std::string>("od");
  const std::string outfile = clp.getValue<std::string>("o");
  std::string earliestTime = clp.getValue<std::string>("et");
  std::string latestTime = clp.getValue<std::string>("lt");

  const int earliestDep = secondsSinceMonMidnight(earliestDay, earliestTime);
  const int latestDep = secondsSinceMonMidnight(latestDay, latestTime);
  if (latestDep < earliestDep)
    throw std::invalid_argument("latest departure before earliest departure");

  std::cout << "Reading the vertex sets..." << std::flush;
  std::ifstream in(vertexSetFile, std::ios::binary);
  if (!in.good())
    throw std::invalid_argument("file not found -- '" + vertexSetFile);

  int numZones;
  bio::read(in, numZones);
  if (numZones < 0)
    throw std::invalid_argument("vertex set file corrupt");

  std::unordered_map<int, int> origToSeqZoneIds;
  for (int seqId = 1; seqId <= numZones; ++seqId) {
    int origId;
    bio::read(in, origId);
    origToSeqZoneIds[origId] = seqId;
  }

  std::vector<std::vector<int>> zoneVertexSets(numZones);
  for (auto& vertexSet : zoneVertexSets)
    bio::read(in, vertexSet);

  in.close();
  std::cout << " done." << std::endl;

  std::cout << "Computing the OD-pairs..." << std::flush;
  std::ofstream out(outfile + ".csv");
  if (!out.good())
    throw std::invalid_argument("file cannot be opened -- '" + outfile + ".csv'");
  out << "# Input graph: " << graphFile << std::endl;
  out << "# Methodology: mobiTopp" << std::endl;
  out << "origin,destination,origin_zone,destination_zone,departure" << std::endl;

  std::default_random_engine rand(seed);
  int householdId, personId, travelTime, activity, durOfActivity;
  char *dayOfWeek = nullptr, *dep = nullptr, *mode = nullptr;
  char *arr, *originZone, *destinationZone, *length;
  io::CSVReader<12, io::trim_chars<>, io::no_quote_escape<';'>> od(odFile);
  od.set_header(
      "household_id", "person_id", "day_of_week", "dep_time", "arr_time", "travel_time",
      "origin_zone", "destination_zone", "length", "mode", "activity", "duration_of_activity");
  while (od.read_row(householdId, personId, dayOfWeek, dep, arr, travelTime,
      originZone, destinationZone, length, mode, activity, durOfActivity)) {
    // Does the current OD-pair lie in the selected analysis period?
    const int departure = secondsSinceMonMidnight(dayOfWeek, dep);
    if (departure < earliestDep || latestDep <= departure)
      continue;

    // Does the current OD-pair use the selected mode of transportation?
    if (std::strcmp(mode, selectedMode.c_str()))
      continue;

    if (originZone[0] != 'Z')
      throw std::domain_error("invalid zone ID -- '" + std::string(originZone) + "'");
    if (destinationZone[0] != 'Z')
      throw std::domain_error("invalid zone ID -- '" + std::string(destinationZone) + "'");

    int oZoneId = origToSeqZoneIds[lexicalCast<int>(++originZone)];
    int dZoneId = origToSeqZoneIds[lexicalCast<int>(++destinationZone)];
    if (oZoneId-- == 0)
      throw std::domain_error("zone not found -- '" + std::string(originZone) + "'");
    if (dZoneId-- == 0)
      throw std::domain_error("zone not found -- '" + std::string(destinationZone) + "'");

    const std::vector<int>& originVertexSet = zoneVertexSets[oZoneId];
    const std::vector<int>& destinationVertexSet = zoneVertexSets[dZoneId];
    if (originVertexSet.size() == 0 || destinationVertexSet.size() == 0)
      continue;

    std::uniform_int_distribution<> originDist(0, originVertexSet.size() - 1);
    std::uniform_int_distribution<> destinationDist(0, destinationVertexSet.size() - 1);
    const int o = originVertexSet[originDist(rand)];
    const int d = destinationVertexSet[destinationDist(rand)];
    out << o << ',' << d << ',' << oZoneId << ',' << dZoneId << ',' << departure << std::endl;
  }
  std::cout << " done." << std::endl;
}

int main(int argc, char* argv[]) {
  try {
    CommandLineParser clp(argc, argv);
    if (clp.isSet("help"))
      printUsage();
    else if (clp.isSet("g") && clp.isSet("v") && clp.isSet("o"))
      computeVertexSets(clp);
    else if (clp.isSet("g") && clp.isSet("vs") && clp.isSet("od") && clp.isSet("o"))
      extractODPairs(clp);
    else
      throw std::invalid_argument("invalid options");
  } catch (std::exception& e) {
    std::cerr << argv[0] << ": " << e.what() << std::endl;
    std::cerr << "Try '" << argv[0] <<" -help' for more information." << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
