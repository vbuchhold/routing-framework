#include <cassert>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <csv.h>
#include <routingkit/contraction_hierarchy.h>
#include <routingkit/customizable_contraction_hierarchy.h>
#include <routingkit/nested_dissection.h>

#include "Algorithms/CCH/EliminationTreeQuery.h"
#include "Algorithms/CH/CHConversion.h"
#include "Algorithms/CH/CHPreprocessing.h"
#include "Algorithms/CH/CHQuery.h"
#include "Algorithms/CH/ContractionHierarchy.h"
#include "Algorithms/Dijkstra/BiDijkstra.h"
#include "Algorithms/Dijkstra/Dijkstra.h"
#include "DataStructures/Graph/Attributes/EdgeIdAttribute.h"
#include "DataStructures/Graph/Attributes/LatLngAttribute.h"
#include "DataStructures/Graph/Attributes/LengthAttribute.h"
#include "DataStructures/Graph/Attributes/TravelTimeAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "Tools/CommandLine/CommandLineParser.h"
#include "Tools/BinaryIO.h"
#include "Tools/Constants.h"
#include "Tools/Timer.h"

void printUsage() {
  std::cout <<
      "Usage: RunP2PAlgo -a CH -g <file> [-ord <file>] -o <file>\n"
      "       RunP2PAlgo -a CCH [-b <balance>] -g <file> -o <file>\n"
      "       RunP2PAlgo -a Dij -g <file> -od <file> -o <file>\n"
      "       RunP2PAlgo -a Bi-Dij -g <file> -od <file> -o <file>\n"
      "       RunP2PAlgo -a CH [-s] -ch <file> -od <file> -o <file>\n"
      "       RunP2PAlgo -a CCH-Dij [-s] -g <file> -ord <file> -od <file> -o <file>\n"
      "       RunP2PAlgo -a CCH-tree -g <file> -ord <file> -od <file> -o <file>\n"
      "This program runs the preprocessing or query phase of various point-to-point\n"
      "shortest-path algorithms, such as Dijkstra, bidirectional search, CH or CCH.\n"
      "  -l                use physical length as metric (default: travel time)\n"
      "  -s                do not use the stall-on-demand technique\n"
      "  -a <algo>         algorithm to be run\n"
      "  -b <balance>      balance parameter in % for nested dissection (default: 30)\n"
      "  -g <file>         input graph in binary format\n"
      "  -ord <file>       order in which vertices are contracted\n"
      "  -ch <file>        metric-dependent CH\n"
      "  -od <file>        file containing OD-pairs (queries)\n"
      "  -o <file>         place output in <file>\n"
      "  -help             display this help and exit\n";
}

// Some helper aliases.
using VertexAttributes = VertexAttrs<LatLngAttribute>;
using EdgeAttributes = EdgeAttrs<EdgeIdAttribute, LengthAttribute, TravelTimeAttribute>;
using InputGraph = StaticGraph<VertexAttributes, EdgeAttributes>;
using CHGraph = StaticGraph<VertexAttrs<>, EdgeAttrs<EdgeIdAttribute, TravelTimeAttribute>>;
using CH = ContractionHierarchy<CHGraph, TravelTimeAttribute>;
#ifdef P2P_OUTPUT_PHYSICAL_PATH_LENGTHS
using LabelSet = BasicLabelSet<0, ParentInfo::FULL_PARENT_INFO>;
#else
using LabelSet = BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>;
#endif

// The query algorithms.
using Dij = StandardDijkstra<InputGraph, TravelTimeAttribute, LabelSet>;
using BiDij = BiDijkstra<Dij>;
template <bool useStalling>
using CHSearch = StandardCHQuery<CH, LabelSet, useStalling>;
using CCHTree = EliminationTreeQuery<CH, LabelSet>;

#ifdef P2P_OUTPUT_PHYSICAL_PATH_LENGTHS
InputGraph inputGraph;
#endif

// Writes the header line of the output CSV file.
template <typename AlgoT>
inline void writeHeaderLine(std::ofstream& out, AlgoT&) {
  out << "query_time,distance";
#ifdef P2P_OUTPUT_PHYSICAL_PATH_LENGTHS
  out << ",length";
#endif
  out << '\n';
}

// Writes a record line of the output CSV file, containing statistics about a single query.
template <typename AlgoT>
inline void writeRecordLine(std::ofstream& out, AlgoT& algo, const int, const int elapsed) {
  out << elapsed << ',' << algo.getDistance();
#ifdef P2P_OUTPUT_PHYSICAL_PATH_LENGTHS
  int len = 0;
  for (const auto e : algo.getEdgePath())
    len += inputGraph.length(e);
  out << ',' << len;
#endif
  out << '\n';
}

template <>
inline void writeRecordLine(std::ofstream& out, Dij& algo, const int dst, const int elapsed) {
  out << elapsed << ',' << algo.getDistance(dst);
#ifdef P2P_OUTPUT_PHYSICAL_PATH_LENGTHS
  int len = 0;
  for (const auto e : algo.getReverseEdgePath(dst))
    len += inputGraph.length(e);
  out << ',' << len;
#endif
  out << '\n';
}

// Runs the queries in the specified OD-file using the specified P2P algorithm.
template <typename AlgoT, typename T>
inline void runQueries(AlgoT& algo, const std::string& od, std::ofstream& outfile, T translate) {
  Timer timer;
  int origin, destination, rank;
  using TrimPolicy = io::trim_chars<>;
  using QuotePolicy = io::no_quote_escape<','>;
  using OverflowPolicy = io::throw_on_overflow;
  using CommentPolicy = io::single_line_comment<'#'>;
  io::CSVReader<3, TrimPolicy, QuotePolicy, OverflowPolicy, CommentPolicy> odFile(od);
  const io::ignore_column ignore = io::ignore_extra_column | io::ignore_missing_column;
  odFile.read_header(ignore, "origin", "destination", "dijkstra_rank");
  const bool hasRanks = odFile.has_column("dijkstra_rank");
  if (hasRanks) outfile << "dijkstra_rank,";
  writeHeaderLine(outfile, algo);
  while (odFile.read_row(origin, destination, rank)) {
    origin = translate(origin);
    destination = translate(destination);
    timer.restart();
    algo.run(origin, destination);
    const int elapsed = timer.elapsed<std::chrono::microseconds>();
    if (hasRanks) outfile << rank << ',';
    writeRecordLine(outfile, algo, destination, elapsed);
  }
}

// Invoked when the user wants to run the query phase of a P2P algorithm.
inline void runQueries(const CommandLineParser& clp) {
  const bool useLengths = clp.isSet("l");
  const bool noStalling = clp.isSet("s");
  const std::string algorithmName = clp.getValue<std::string>("a");
  const std::string graphFilename = clp.getValue<std::string>("g");
  const std::string orderFilename = clp.getValue<std::string>("ord");
  const std::string chFilename = clp.getValue<std::string>("ch");
  const std::string odFilename = clp.getValue<std::string>("od");
  const std::string outfilename = clp.getValue<std::string>("o");

#ifdef P2P_OUTPUT_PHYSICAL_PATH_LENGTHS
  // Read the input graph.
  std::ifstream graphFile(graphFilename, std::ios::binary);
  if (!graphFile.good())
    throw std::invalid_argument("file not found -- '" + graphFilename + "'");
  inputGraph.readFrom(graphFile);
  graphFile.close();
  int id = 0;
  FORALL_EDGES(inputGraph, v)
    inputGraph.edgeId(v) = id++;
#endif

  // Open the output CSV file.
  std::ofstream outfile(outfilename + ".csv");
  if (!outfile.good())
    throw std::invalid_argument("file cannot be opened -- '" + outfilename + ".csv'");

  if (algorithmName == "Dij") {
    // Run the query phase of Dijkstra's algorithm.
    std::ifstream graphFile(graphFilename, std::ios::binary);
    if (!graphFile.good())
      throw std::invalid_argument("file not found -- '" + graphFilename + "'");
    InputGraph graph(graphFile);
    graphFile.close();

    if (useLengths)
      FORALL_EDGES(graph, e)
        graph.travelTime(e) = graph.length(e);

    outfile << "# graph: " << graphFilename << '\n';
    outfile << "# OD-pairs: " << odFilename << '\n';

    Dij algo(graph);
    runQueries(algo, odFilename, outfile, [](const int v) { return v; });
  } else if (algorithmName == "Bi-Dij") {
    // Run the query phase of bidirectional search.
    std::ifstream graphFile(graphFilename, std::ios::binary);
    if (!graphFile.good())
      throw std::invalid_argument("file not found -- '" + graphFilename + "'");
    InputGraph graph(graphFile);
    graphFile.close();

    if (useLengths)
      FORALL_EDGES(graph, e)
        graph.travelTime(e) = graph.length(e);

    outfile << "# graph: " << graphFilename << '\n';
    outfile << "# OD-pairs: " << odFilename << '\n';

    InputGraph reverse = graph.getReverseGraph();
    BiDij algo(graph, reverse);
    runQueries(algo, odFilename, outfile, [](const int v) { return v; });
  } else if (algorithmName == "CH") {
    // Run the query phase of CH.
    std::ifstream chFile(chFilename, std::ios::binary);
    if (!chFile.good())
      throw std::invalid_argument("file not found -- '" + chFilename + "'");
    CH ch(chFile);
    chFile.close();

    outfile << "# CH: " << chFilename << '\n';
    outfile << "# OD-pairs: " << odFilename << '\n';

    if (noStalling) {
      CHSearch<false> algo(ch);
      runQueries(algo, odFilename, outfile, [&](const int v) { return ch.rank(v); });
    } else {
      CHSearch<true> algo(ch);
      runQueries(algo, odFilename, outfile, [&](const int v) { return ch.rank(v); });
    }
  } else if (algorithmName == "CCH-Dij") {
    // Run the Dijkstra-based query phase of CCH.
    std::ifstream graphFile(graphFilename, std::ios::binary);
    if (!graphFile.good())
      throw std::invalid_argument("file not found -- '" + graphFilename + "'");
    InputGraph graph(graphFile);
    graphFile.close();

    std::ifstream orderFile(orderFilename, std::ios::binary);
    if (!orderFile.good())
      throw std::invalid_argument("file not found -- '" + orderFilename + "'");
    std::vector<unsigned int> order;
    bio::read(orderFile, order);
    orderFile.close();

    std::vector<unsigned int> tails(graph.numEdges());
    std::vector<unsigned int> heads(graph.numEdges());
    FORALL_VERTICES(graph, u)
      FORALL_INCIDENT_EDGES(graph, u, e) {
        tails[e] = u;
        heads[e] = graph.edgeHead(e);
      }

    using RoutingKitCCH = RoutingKit::CustomizableContractionHierarchy;
    using RoutingKitCCHMetric = RoutingKit::CustomizableContractionHierarchyMetric;
    using RoutingKitCH = RoutingKit::ContractionHierarchy;
    RoutingKitCCH cch(order, tails, heads);
    const int* const weights = useLengths ? &graph.length(0) : &graph.travelTime(0);
    RoutingKitCCHMetric metric(cch, reinterpret_cast<const unsigned int*>(weights));
    RoutingKitCH ch = metric.build_contraction_hierarchy_using_perfect_witness_search();
    CH perfectCH = convert<CH>(ch, graph.numEdges());

    outfile << "# graph: " << graphFilename << '\n';
    outfile << "# order: " << orderFilename << '\n';
    outfile << "# OD-pairs: " << odFilename << '\n';

    if (noStalling) {
      CHSearch<false> algo(perfectCH);
      runQueries(algo, odFilename, outfile, [&](const int v) { return perfectCH.rank(v); });
    } else {
      CHSearch<true> algo(perfectCH);
      runQueries(algo, odFilename, outfile, [&](const int v) { return perfectCH.rank(v); });
    }
  } else if (algorithmName == "CCH-tree") {
    // Run the elimination-tree-based query phase of CCH.
    std::ifstream graphFile(graphFilename, std::ios::binary);
    if (!graphFile.good())
      throw std::invalid_argument("file not found -- '" + graphFilename + "'");
    InputGraph graph(graphFile);
    graphFile.close();

    std::ifstream orderFile(orderFilename, std::ios::binary);
    if (!orderFile.good())
      throw std::invalid_argument("file not found -- '" + orderFilename + "'");
    std::vector<unsigned int> order;
    bio::read(orderFile, order);
    orderFile.close();

    std::vector<unsigned int> tails(graph.numEdges());
    std::vector<unsigned int> heads(graph.numEdges());
    FORALL_VERTICES(graph, u)
      FORALL_INCIDENT_EDGES(graph, u, e) {
        tails[e] = u;
        heads[e] = graph.edgeHead(e);
      }

    using RoutingKitCCH = RoutingKit::CustomizableContractionHierarchy;
    using RoutingKitCCHMetric = RoutingKit::CustomizableContractionHierarchyMetric;
    using RoutingKitCH = RoutingKit::ContractionHierarchy;
    RoutingKitCCH cch(order, tails, heads);
    std::vector<int> tree(cch.elimination_tree_parent.begin(), cch.elimination_tree_parent.end());
    tree.back() = INVALID_VERTEX;
    const int* const weights = useLengths ? &graph.length(0) : &graph.travelTime(0);
    RoutingKitCCHMetric metric(cch, reinterpret_cast<const unsigned int*>(weights));
    RoutingKitCH ch = metric.build_contraction_hierarchy_using_perfect_witness_search();
    CH perfectCH = convert<CH>(ch, graph.numEdges());

    outfile << "# graph: " << graphFilename << '\n';
    outfile << "# order: " << orderFilename << '\n';
    outfile << "# OD-pairs: " << odFilename << '\n';

    CCHTree algo(perfectCH, tree);
    runQueries(algo, odFilename, outfile, [&](const int v) { return perfectCH.rank(v); });
  } else {
    throw std::invalid_argument("invalid P2P algorithm -- '" + algorithmName + "'");
  }
}

// Invoked when the user wants to run the preprocessing phase of a P2P algorithm.
inline void runPreprocessing(const CommandLineParser& clp) {
  const bool useLengths = clp.isSet("l");
  const int imbalance = clp.getValue<int>("i", 30);
  const std::string algorithmName = clp.getValue<std::string>("a");
  const std::string graphFilename = clp.getValue<std::string>("g");
  const std::string orderFilename = clp.getValue<std::string>("ord");
  std::string outfilename = clp.getValue<std::string>("o");

  // Read the input graph.
  std::ifstream graphFile(graphFilename, std::ios::binary);
  if (!graphFile.good())
    throw std::invalid_argument("file not found -- '" + graphFilename + "'");
  InputGraph graph(graphFile);
  graphFile.close();

  if (algorithmName == "CH") {
    // Run the preprocessing phase of CH.
    if (useLengths)
      FORALL_EDGES(graph, e)
        graph.travelTime(e) = graph.length(e);

    outfilename += (useLengths ? ".dist" : ".time") + std::string(".ch.bin");
    std::ofstream outfile(outfilename, std::ios::binary);
    if (!outfile.good())
      throw std::invalid_argument("file cannot be opened -- '" + outfilename);

    CHPreprocessing<InputGraph, CHGraph, TravelTimeAttribute> algo(graph);
    algo.run().writeTo(outfile);
  } else if (algorithmName == "CCH") {
    // Run the preprocessing phase of CCH.
    if (imbalance < 0)
      throw std::invalid_argument("invalid imbalance -- '" + std::to_string(imbalance) + "'");

    std::vector<float> lats(graph.numVertices());
    std::vector<float> lngs(graph.numVertices());
    std::vector<unsigned int> tails(graph.numEdges());
    std::vector<unsigned int> heads(graph.numEdges());
    FORALL_VERTICES(graph, u) {
      lats[u] = graph.latLng(u).latInDeg();
      lngs[u] = graph.latLng(u).lngInDeg();
      FORALL_INCIDENT_EDGES(graph, u, e) {
        tails[e] = u;
        heads[e] = graph.edgeHead(e);
      }
    }

    const RoutingKit::GraphFragment fragment =
        RoutingKit::make_graph_fragment(graph.numVertices(), tails, heads);

    auto computeSeparator = [&lats, &lngs, imbalance](const RoutingKit::GraphFragment& fragment) {
      const RoutingKit::CutSide side = inertial_flow(fragment, imbalance, lats, lngs);
      return derive_separator_from_cut(fragment, side.is_node_on_side);
    };

    outfilename += ".nd" + std::to_string(imbalance) + ".ord.bin";
    std::ofstream outfile(outfilename, std::ios::binary);
    if (!outfile.good())
      throw std::invalid_argument("file cannot be opened -- '" + outfilename);
    bio::write(outfile, compute_nested_node_dissection_order(fragment, computeSeparator));
  } else {
    throw std::invalid_argument("invalid P2P algorithm -- '" + algorithmName + "'");
  }
}

int main(int argc, char* argv[]) {
  try {
    CommandLineParser clp(argc, argv);
    if (clp.isSet("help"))
      printUsage();
    else if (clp.isSet("od"))
      runQueries(clp);
    else
      runPreprocessing(clp);
  } catch (std::exception& e) {
    std::cerr << argv[0] << ": " << e.what() << std::endl;
    std::cerr << "Try '" << argv[0] <<" -help' for more information." << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
