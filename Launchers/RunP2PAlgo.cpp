#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <csv.h>
#include <routingkit/nested_dissection.h>

#include "Algorithms/CCH/CCH.h"
#include "Algorithms/CCH/CCHMetric.h"
#include "Algorithms/CCH/EliminationTreeQuery.h"
#include "Algorithms/CH/CH.h"
#include "Algorithms/CH/CHQuery.h"
#include "Algorithms/Dijkstra/BiDijkstra.h"
#include "Algorithms/Dijkstra/Dijkstra.h"
#include "DataStructures/Graph/Attributes/LatLngAttribute.h"
#include "DataStructures/Graph/Attributes/LengthAttribute.h"
#include "DataStructures/Graph/Attributes/TravelTimeAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "DataStructures/Partitioning/SeparatorDecomposition.h"
#include "Tools/CommandLine/CommandLineParser.h"
#include "Tools/Timer.h"

void printUsage() {
  std::cout <<
      "Usage: RunP2PAlgo -a CH -g <file> -o <file>\n"
      "       RunP2PAlgo -a CCH [-b <balance>] -g <file> -o <file>\n"
      "       RunP2PAlgo -a CCH-custom -n <num> -g <file> -sep <file> -o <file>\n"
      "       RunP2PAlgo -a Dij -g <file> -od <file> -o <file>\n"
      "       RunP2PAlgo -a Bi-Dij -g <file> -od <file> -o <file>\n"
      "       RunP2PAlgo -a CH [-s] -ch <file> -od <file> -o <file>\n"
      "       RunP2PAlgo -a CCH-Dij [-s] -g <file> -sep <file> -od <file> -o <file>\n"
      "       RunP2PAlgo -a CCH-tree -g <file> -sep <file> -od <file> -o <file>\n"
      "This program runs the preprocessing or query phase of various point-to-point\n"
      "shortest-path algorithms, such as Dijkstra, bidirectional search, CH or CCH.\n"
      "  -l                use physical length as metric (default: travel time)\n"
      "  -s                do not use the stall-on-demand technique\n"
      "  -a <algo>         algorithm to be run\n"
      "  -b <balance>      balance parameter in % for nested dissection (default: 30)\n"
      "  -n <num>          run the customization <num> times (defaults to 100)\n"
      "  -g <file>         input graph in binary format\n"
      "  -ch <file>        weighted CH\n"
      "  -sep <file>       separator decomposition of the input graph\n"
      "  -od <file>        file containing OD-pairs (queries)\n"
      "  -o <file>         place output in <file>\n"
      "  -help             display this help and exit\n";
}

// Some helper aliases.
using VertexAttributes = VertexAttrs<LatLngAttribute>;
using EdgeAttributes = EdgeAttrs<LengthAttribute, TravelTimeAttribute>;
using InputGraph = StaticGraph<VertexAttributes, EdgeAttributes>;
using LabelSet = BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>;

// The query algorithms.
using Dij = StandardDijkstra<InputGraph, TravelTimeAttribute, LabelSet>;
using BiDij = BiDijkstra<Dij>;
template <bool useStalling>
using CHSearch = StandardCHQuery<LabelSet, useStalling>;
using CCHTree = EliminationTreeQuery<LabelSet>;

// Writes the header line of the output CSV file.
template <typename AlgoT>
inline void writeHeaderLine(std::ofstream& out, AlgoT&) {
  out << "query_time,distance" << '\n';
}

// Writes a record line of the output CSV file, containing statistics about a single query.
template <typename AlgoT>
inline void writeRecordLine(std::ofstream& out, AlgoT& algo, const int, const int elapsed) {
  out << elapsed << ',' << algo.getDistance() << '\n';
}

template <>
inline void writeRecordLine(std::ofstream& out, Dij& algo, const int dst, const int elapsed) {
  out << elapsed << ',' << algo.getDistance(dst) << '\n';
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
    const auto elapsed = timer.elapsed<std::chrono::microseconds>();
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
  const std::string orderFilename = clp.getValue<std::string>("sep");
  const std::string chFilename = clp.getValue<std::string>("ch");
  const std::string odFilename = clp.getValue<std::string>("od");
  const std::string outfilename = clp.getValue<std::string>("o");

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
    SeparatorDecomposition sepDecomp;
    sepDecomp.readFrom(orderFile);
    orderFile.close();

    CCH cch;
    cch.preprocess(graph, sepDecomp);
    CCHMetric metric(cch, useLengths ? &graph.length(0) : &graph.travelTime(0));
    const auto minWeightedCH = metric.buildMinimumWeightedCH();

    outfile << "# graph: " << graphFilename << '\n';
    outfile << "# order: " << orderFilename << '\n';
    outfile << "# OD-pairs: " << odFilename << '\n';

    if (noStalling) {
      CHSearch<false> algo(minWeightedCH);
      runQueries(algo, odFilename, outfile, [&](const int v) { return minWeightedCH.rank(v); });
    } else {
      CHSearch<true> algo(minWeightedCH);
      runQueries(algo, odFilename, outfile, [&](const int v) { return minWeightedCH.rank(v); });
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
    SeparatorDecomposition sepDecomp;
    sepDecomp.readFrom(orderFile);
    orderFile.close();

    CCH cch;
    cch.preprocess(graph, sepDecomp);
    CCHMetric metric(cch, useLengths ? &graph.length(0) : &graph.travelTime(0));
    const auto minWeightedCH = metric.buildMinimumWeightedCH();

    outfile << "# graph: " << graphFilename << '\n';
    outfile << "# order: " << orderFilename << '\n';
    outfile << "# OD-pairs: " << odFilename << '\n';

    CCHTree algo(minWeightedCH, cch.getEliminationTree());
    runQueries(algo, odFilename, outfile, [&](const int v) { return minWeightedCH.rank(v); });
  } else {
    throw std::invalid_argument("invalid P2P algorithm -- '" + algorithmName + "'");
  }
}

// Invoked when the user wants to run the preprocessing phase of a P2P algorithm.
inline void runPreprocessing(const CommandLineParser& clp) {
  const auto useLengths = clp.isSet("l");
  const auto imbalance = clp.getValue<int>("i", 30);
  const auto numCustomRuns = clp.getValue<int>("n", 100);
  const auto algorithmName = clp.getValue<std::string>("a");
  const auto graphFilename = clp.getValue<std::string>("g");
  const auto orderFilename = clp.getValue<std::string>("sep");
  auto outfilename = clp.getValue<std::string>("o");

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

    CH ch;
    ch.preprocess<TravelTimeAttribute>(graph);
    ch.writeTo(outfile);
  } else if (algorithmName == "CCH") {
    // Run the preprocessing phase of CCH.
    if (imbalance < 0)
      throw std::invalid_argument("invalid imbalance -- '" + std::to_string(imbalance) + "'");

    // Convert the input graph to RoutingKit's graph representation.
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

    // Compute a separator decomposition for the input graph.
    const auto fragment = RoutingKit::make_graph_fragment(graph.numVertices(), tails, heads);
    auto computeSep = [&](const RoutingKit::GraphFragment& fragment) {
      const auto cut = inertial_flow(fragment, imbalance, lats, lngs);
      return derive_separator_from_cut(fragment, cut.is_node_on_side);
    };
    const auto decomp = compute_separator_decomposition(fragment, computeSep);

    // Convert the separator decomposition to our representation.
    SeparatorDecomposition sepDecomp;
    for (const auto& n : decomp.tree) {
      SeparatorDecomposition::Node node;
      node.leftChild = n.left_child;
      node.rightSibling = n.right_sibling;
      node.firstSeparatorVertex = n.first_separator_vertex;
      node.lastSeparatorVertex = n.last_separator_vertex;
      sepDecomp.tree.push_back(node);
    }
    sepDecomp.order.assign(decomp.order.begin(), decomp.order.end());

    outfilename += ".nd" + std::to_string(imbalance) + ".sep.bin";
    std::ofstream outfile(outfilename, std::ios::binary);
    if (!outfile.good())
      throw std::invalid_argument("file cannot be opened -- '" + outfilename);
    sepDecomp.writeTo(outfile);
  } else if (algorithmName == "CCH-custom") {
    // Run the customization phase of CCH.
    std::ifstream orderFile(orderFilename, std::ios::binary);
    if (!orderFile.good())
      throw std::invalid_argument("file not found -- '" + orderFilename + "'");
    SeparatorDecomposition decomp;
    decomp.readFrom(orderFile);
    orderFile.close();

    std::ofstream outfile(outfilename + ".csv");
    if (!outfile.good())
      throw std::invalid_argument("file cannot be opened -- '" + outfilename + ".csv'");
    outfile << "# graph: " << graphFilename << '\n';
    outfile << "# order: " << orderFilename << '\n';
    outfile << "basic_customization,perfect_customization,construction,total_time\n";

    CCH cch;
    cch.preprocess(graph, decomp);

    Timer timer;
    int basicCustom, perfectCustom, construction, total;
    for (auto i = 0; i < numCustomRuns; ++i) {
      {
        CCHMetric metric(cch, &graph.travelTime(0));
        timer.restart();
        metric.customize();
        basicCustom = timer.elapsed<std::chrono::microseconds>();
        timer.restart();
        metric.runPerfectCustomization();
        perfectCustom = timer.elapsed<std::chrono::microseconds>();
      }
      {
        CCHMetric metric(cch, &graph.travelTime(0));
        timer.restart();
        metric.buildMinimumWeightedCH();
        total = timer.elapsed<std::chrono::microseconds>();
      }
      construction = total - basicCustom - perfectCustom;
      outfile << basicCustom << ',' << perfectCustom << ',' << construction << ',' << total << '\n';
    }
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
