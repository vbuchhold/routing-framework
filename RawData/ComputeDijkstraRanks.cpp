#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

#include <csv.h>

#include "DataStructures/Graph/Attributes/LengthAttribute.h"
#include "DataStructures/Graph/Attributes/TravelTimeAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "Experiments/ODPairGenerator.h"
#include "Tools/CommandLine/CommandLineParser.h"
#include "Tools/Constants.h"

void printUsage() {
  std::cout <<
      "Usage: ComputeDijkstraRanks -od <file> -g <file> -o <file>\n"
      "This program computes the Dijkstra ranks for the specified OD-pairs.\n"
      "  -l                compute ranks wrt physical lengths (default: travel times)\n"
      "  -od <file>        input OD-file\n"
      "  -g <file>         input graph in binary format\n"
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
    const std::string graphFilename = clp.getValue<std::string>("g");
    const std::string outfilename = clp.getValue<std::string>("o");

    using Graph = StaticGraph<VertexAttrs<>, EdgeAttrs<LengthAttribute, TravelTimeAttribute>>;
    std::ifstream graphFile(graphFilename, std::ios::binary);
    if (!graphFile.good())
      throw std::invalid_argument("file not found -- '" + graphFilename + "'");
    Graph graph(graphFile);
    graphFile.close();
    std::default_random_engine rand;
    ODPairGenerator<Graph, TravelTimeAttribute> gen(graph, rand);

    if (clp.isSet("l"))
      // Compute Dijkstra ranks wrt physical lengths.
      FORALL_EDGES(graph, e)
        graph.travelTime(e) = graph.length(e);

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

    int origin = INVALID_VERTEX, destination = INVALID_VERTEX, originZone, destinationZone, dep;
    using TrimPolicy = io::trim_chars<>;
    using QuotePolicy = io::no_quote_escape<','>;
    using OverflowPolicy = io::throw_on_overflow;
    using CommentPolicy = io::single_line_comment<'#'>;
    io::CSVReader<5, TrimPolicy, QuotePolicy, OverflowPolicy, CommentPolicy> csv(odFilename);
    csv.read_header(
        io::ignore_extra_column | io::ignore_missing_column,
        "origin", "destination", "origin_zone", "destination_zone", "departure");
    if (csv.has_column("origin_zone") != csv.has_column("destination_zone"))
      throw std::invalid_argument("OD-file corrupt");
    const bool hasZones = csv.has_column("origin_zone");
    const bool hasDep = csv.has_column("departure");

    out << "origin,destination";
    if (hasZones) out << ",origin_zone,destination_zone";
    if (hasDep) out << ",departure";
    out << ",dijkstra_rank\n";

    while (csv.read_row(origin, destination, originZone, destinationZone, dep)) {
      out << origin << ',' << destination;
      if (hasZones) out << ',' << originZone << ',' << destinationZone;
      if (hasDep) out << ',' << dep;
      out << ',' << gen.getDijkstraRankFor({origin, destination}) << '\n';
    }
  } catch (std::exception& e) {
    std::cerr << argv[0] << ": " << e.what() << std::endl;
    std::cerr << "Try '" << argv[0] <<" -help' for more information." << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
