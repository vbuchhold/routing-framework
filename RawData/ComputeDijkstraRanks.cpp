#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <csv.h>

#include "DataStructures/Graph/Attributes/LengthAttribute.h"
#include "DataStructures/Graph/Attributes/SequentialVertexIdAttribute.h"
#include "DataStructures/Graph/Attributes/TravelTimeAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "Experiments/ODPairGenerator.h"
#include "Tools/CommandLine/CommandLineParser.h"
#include "Tools/Constants.h"

void printUsage() {
  std::cout <<
      "Usage: ComputeDijkstraRanks [-l] -od <file> -g <file> -o <file>\n"
      "Compute Dijkstra ranks for given OD-pairs.\n"
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

    // Read the graph from file.
    using VertexAttributes = VertexAttrs<SequentialVertexIdAttribute>;
    using EdgeAttributes = EdgeAttrs<LengthAttribute, TravelTimeAttribute>;
    using Graph = StaticGraph<VertexAttributes, EdgeAttributes>;
    std::ifstream graphFile(graphFilename, std::ios::binary);
    if (!graphFile.good())
      throw std::invalid_argument("file not found -- '" + graphFilename + "'");
    Graph graph(graphFile);
    graphFile.close();

    if (clp.isSet("l"))
      // Compute Dijkstra ranks wrt physical lengths.
      FORALL_EDGES(graph, e)
        graph.travelTime(e) = graph.length(e);

    // Build a map to translate an original into a local vertex identifier.
    std::vector<int> origToLocalId(graph.numVertices(), INVALID_VERTEX);
    if (graph.numVertices() > 0 && graph.sequentialVertexId(0) == INVALID_VERTEX) {
      FORALL_VERTICES(graph, v)
        origToLocalId[v] = v;
    } else {
      FORALL_VERTICES(graph, v) {
        const int origId = graph.sequentialVertexId(v);
        if (origToLocalId.size() < origId + 1)
          origToLocalId.resize(origId + 1, INVALID_VERTEX);
        origToLocalId[origId] = v;
      }
    }

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
        "origin", "destination", "origin_zone", "destination_zone", "departure", "distance");
    if (!csv.has_column("origin") || !csv.has_column("destination"))
      throw std::invalid_argument("OD-file corrupt");
    if (csv.has_column("origin_zone") != csv.has_column("destination_zone"))
      throw std::invalid_argument("OD-file corrupt");
    const bool hasZones = csv.has_column("origin_zone");
    const bool hasDep = csv.has_column("departure");
    const bool hasDist = csv.has_column("distance");

    // Write the header line to the output OD-file.
    out << "origin,destination";
    if (hasZones) out << ",origin_zone,destination_zone";
    if (hasDep) out << ",departure";
    out << ",dijkstra_rank";
    if (hasDist) out << ",distance";
    out << '\n';

    // Compute the Dijkstra rank for each OD-pair.
    ODPairGenerator<Graph, TravelTimeAttribute> gen(graph);
    int origin, destination, originZone, destinationZone, dep, dist;
    while (csv.read_row(origin, destination, originZone, destinationZone, dep, dist)) {
      out << origin << ',' << destination;
      if (hasZones) out << ',' << originZone << ',' << destinationZone;
      if (hasDep) out << ',' << dep;
      out << ',' << gen.getDijkstraRankFor({origToLocalId[origin], origToLocalId[destination]});
      if (hasDist) out << ',' << dist;
      out << '\n';
    }
  } catch (std::exception& e) {
    std::cerr << argv[0] << ": " << e.what() << std::endl;
    std::cerr << "Try '" << argv[0] <<" -help' for more information." << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
