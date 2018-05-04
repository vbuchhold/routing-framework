#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include <csv.h>

#include "Algorithms/CH/CHQuery.h"
#include "Algorithms/CH/ContractionHierarchy.h"
#include "DataStructures/Graph/Attributes/EdgeIdAttribute.h"
#include "DataStructures/Graph/Attributes/TravelTimeAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "Tools/CommandLine/CommandLineParser.h"

void printUsage() {
  std::cout <<
      "Usage: ComputeDistances -od <file> -ch <file> -o <file>\n"
      "Compute distances for given OD-pairs.\n"
      "  -od <file>        input OD-file\n"
      "  -ch <file>        input CH in binary format\n"
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

    const std::string odFilename = clp.getValue<std::string>("od");
    const std::string chFilename = clp.getValue<std::string>("ch");
    const std::string outfilename = clp.getValue<std::string>("o");

    // Read the CH from file.
    using CHGraph = StaticGraph<VertexAttrs<>, EdgeAttrs<EdgeIdAttribute, TravelTimeAttribute>>;
    using CH = ContractionHierarchy<CHGraph, TravelTimeAttribute>;
    std::ifstream chFile(chFilename, std::ios::binary);
    if (!chFile.good())
      throw std::invalid_argument("file not found -- '" + chFilename + "'");
    CH ch(chFile);
    chFile.close();

    // Copy the comment lines at the beginning of the input OD-file into the output OD-file.
    std::ifstream in(odFilename);
    if (!in.good())
      throw std::invalid_argument("file not found -- '" + odFilename + "'");
    std::ofstream out(outfilename + ".csv");
    if (!out.good())
      throw std::invalid_argument("file cannot be opened -- '" + outfilename + ".csv'");
    while (in.peek() == '#') {
      std::string comment;
      getline(in, comment);
      out << comment << '\n';
    }
    in.close();

    // Determine which optional columns (e.g., zones) are present in the input OD-file.
    using TrimPolicy = io::trim_chars<>;
    using QuotePolicy = io::no_quote_escape<','>;
    using OverflowPolicy = io::throw_on_overflow;
    using CommentPolicy = io::single_line_comment<'#'>;
    io::CSVReader<6, TrimPolicy, QuotePolicy, OverflowPolicy, CommentPolicy> csv(odFilename);
    csv.read_header(
        io::ignore_extra_column | io::ignore_missing_column,
        "origin", "destination", "origin_zone", "destination_zone", "departure", "dijkstra_rank");
    if (!csv.has_column("origin") || !csv.has_column("destination"))
      throw std::invalid_argument("OD-file corrupt");
    if (csv.has_column("origin_zone") != csv.has_column("destination_zone"))
      throw std::invalid_argument("OD-file corrupt");
    const bool hasZones = csv.has_column("origin_zone");
    const bool hasDep = csv.has_column("departure");
    const bool hasRank = csv.has_column("dijkstra_rank");

    // Write the header line to the output OD-file.
    out << "origin,destination";
    if (hasZones) out << ",origin_zone,destination_zone";
    if (hasDep) out << ",departure";
    if (hasRank) out << ",dijkstra_rank";
    out << ",distance";
    out << '\n';

    // Compute the distance for each OD-pair.
    StandardCHQuery<CH, BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>> chSearch(ch);
    int origin, destination, originZone, destinationZone, dep, rank;
    while (csv.read_row(origin, destination, originZone, destinationZone, dep, rank)) {
      chSearch.run(ch.rank(origin), ch.rank(destination));
      out << origin << ',' << destination;
      if (hasZones) out << ',' << originZone << ',' << destinationZone;
      if (hasDep) out << ',' << dep;
      if (hasRank) out << ',' << rank;
      out << ',' << chSearch.getDistance();
      out << '\n';
    }
  } catch (std::exception& e) {
    std::cerr << argv[0] << ": " << e.what() << std::endl;
    std::cerr << "Try '" << argv[0] <<" -help' for more information." << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
