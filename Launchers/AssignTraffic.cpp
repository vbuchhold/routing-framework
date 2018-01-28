#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include <routingkit/customizable_contraction_hierarchy.h>
#include <routingkit/nested_dissection.h>

#include "Algorithms/TrafficAssignment/Adapters/BiDijkstraAdapter.h"
#include "Algorithms/TrafficAssignment/Adapters/CCHAdapter.h"
#include "Algorithms/TrafficAssignment/Adapters/CHAdapter.h"
#include "Algorithms/TrafficAssignment/Adapters/DijkstraAdapter.h"
#include "Algorithms/TrafficAssignment/ObjectiveFunctions/SystemOptimum.h"
#include "Algorithms/TrafficAssignment/ObjectiveFunctions/UserEquilibrium.h"
#include "Algorithms/TrafficAssignment/TravelCostFunctions/BprFunction.h"
#include "Algorithms/TrafficAssignment/TravelCostFunctions/DavidsonFunction.h"
#include "Algorithms/TrafficAssignment/TravelCostFunctions/InverseFunction.h"
#include "Algorithms/TrafficAssignment/TravelCostFunctions/ModifiedDavidsonFunction.h"
#include "Algorithms/TrafficAssignment/FrankWolfeAssignment.h"
#include "DataStructures/Graph/Attributes/CapacityAttribute.h"
#include "DataStructures/Graph/Attributes/EdgeIdAttribute.h"
#include "DataStructures/Graph/Attributes/LatLngAttribute.h"
#include "DataStructures/Graph/Attributes/LengthAttribute.h"
#include "DataStructures/Graph/Attributes/TravelCostAttribute.h"
#include "DataStructures/Graph/Attributes/TravelTimeAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Utilities/OriginDestination.h"
#include "DataStructures/Utilities/UnionFind.h"
#include "Tools/CommandLine/CommandLineParser.h"
#include "Tools/BinaryIO.h"

void printUsage() {
  std::cout <<
      "Usage: AssignTraffic [-f <func>] [-a <algo>] -i <file> -od <file> [-o <file>]\n"
      "This program assigns OD-pairs onto a network using the Frank-Wolfe method. It\n"
      "supports different objectives, travel cost functions and shortest-path algos.\n"
      "  -so               find the system optimum (default: user equilibrium)\n"
      "  -v                display informative messages\n"
      "  -p <hrs>          the period of analysis in hours (default: 1)\n"
      "  -n <num>          the number of iterations (0 means use stopping criterion)\n"
      "  -f <func>         the travel cost function\n"
      "                      possible values:\n"
      "                        bpr davidson modified_davidson (default) inverse\n"
      "  -a <algo>         the shortest-path algorithm\n"
      "                      possible values: dijkstra (default) bidijkstra ch cch\n"
      "  -ord <order>      the order of the OD-pairs\n"
      "                      possible values: random input (default) sorted\n"
      "  -s <seed>         start the random number generator with <seed>\n"
      "  -U <num>          the maximum height of a subtree (used for ordering OD-pairs)\n"
      "  -si <intervals>   a blank-separated list of sampling intervals\n"
      "  -i <file>         the input graph in binary format\n"
      "  -od <file>        the OD-pairs to be assigned\n"
      "  -o <file>         the output CSV file without file extension\n"
      "  -fp <file>        output the flow pattern after each iteration in <file>\n"
      "  -help             display this help and exit\n";
}

// Assigns origin and destination zones to OD-pairs based on a partition of the elimination tree.
template <typename GraphT>
inline void assignZonesToODPairs(
    std::vector<ClusteredOriginDestination>& odPairs, const GraphT& inGraph, const int maxHeight) {
  // Convert the input graph to RoutingKit's graph representation.
  const int numVertices = inGraph.numVertices();
  std::vector<float> lats(numVertices);
  std::vector<float> lngs(numVertices);
  std::vector<unsigned int> tails(inGraph.numEdges());
  std::vector<unsigned int> heads(inGraph.numEdges());
  FORALL_VERTICES(inGraph, u) {
    lats[u] = inGraph.latLng(u).latInDeg();
    lngs[u] = inGraph.latLng(u).lngInDeg();
    FORALL_INCIDENT_EDGES(inGraph, u, e) {
      tails[e] = u;
      heads[e] = inGraph.edgeHead(e);
    }
  }

  // Compute the contraction order and the corresponding elimination tree.
  const auto graph = RoutingKit::make_graph_fragment(numVertices, tails, heads);
  auto computeSep = [&lats, &lngs](const RoutingKit::GraphFragment& graph) {
    return derive_separator_from_cut(graph, inertial_flow(graph, 30, lats, lngs).is_node_on_side);
  };
  const auto order = compute_nested_node_dissection_order(graph, computeSep);
  const auto cch = RoutingKit::CustomizableContractionHierarchy(order, tails, heads);
  std::vector<int> tree(cch.elimination_tree_parent.begin(), cch.elimination_tree_parent.end());

  // Decompose the elimination tree into as few cells with bounded diameter as possible.
  UnionFind partition(numVertices);
  std::vector<int> height(numVertices, 0); // height[v] is the height of the subtree rooted at v.
  for (int i = 0; i != numVertices - 1; ++i)
    if (height[i] < maxHeight) {
      partition.unite(i, tree[i]);
      height[tree[i]] = std::max(height[tree[i]], height[i] + 1);
    }

  // Assign origin and destination zones to OD-pairs.
  std::vector<int> cellId(numVertices);
  for (int i = 0; i != numVertices; ++i)
    cellId[order[i]] = partition.find(i);
  for (auto& od : odPairs) {
    od.originZone = cellId[od.origin];
    od.destinationZone = cellId[od.destination];
  }
}

// Assigns all OD-flows onto the input graph.
template <typename FrankWolfeAssignmentT>
void assignTraffic(const CommandLineParser& clp) {
  const std::string infile = clp.getValue<std::string>("i");
  const std::string odfile = clp.getValue<std::string>("od");
  const std::string csvfile = clp.getValue<std::string>("o");
  const std::string outfile = clp.getValue<std::string>("fp");
  const std::string ord = clp.getValue<std::string>("ord", "input");
  const int maxHeight = clp.getValue<int>("U", 32);
  const double period = clp.getValue<double>("p", 1);

  std::ifstream in(infile, std::ios::binary);
  if (!in.good())
    throw std::invalid_argument("file not found -- '" + infile + "'");
  typename FrankWolfeAssignmentT::InputGraph graph(in);
  in.close();

  int id = 0;
  FORALL_EDGES(graph, e) {
    graph.capacity(e) = std::max(std::round(period * graph.capacity(e)), 1.0);
    graph.edgeId(e) = id++;
  }

  std::vector<ClusteredOriginDestination> odPairs = importClusteredODPairsFrom(odfile);
  if (ord == "random") {
    std::default_random_engine rand(clp.getValue<int>("s", 19900325));
    std::shuffle(odPairs.begin(), odPairs.end(), rand);
  } else if (ord == "sorted") {
    assignZonesToODPairs(odPairs, graph, maxHeight);
    std::sort(odPairs.begin(), odPairs.end());
  } else if (ord != "input") {
    throw std::invalid_argument("invalid order -- '" + ord + "'");
  }

  const int numIterations = clp.getValue<int>("n");
  if (numIterations < 0) {
    const std::string msg("negative number of iterations");
    throw std::invalid_argument(msg + " -- " + std::to_string(numIterations));
  }

  const auto intervals = clp.getValues<int>("si");
  if (intervals[0] <= 0) {
    const std::string msg("non-positive sampling interval");
    throw std::invalid_argument(msg + " -- " + std::to_string(intervals[0]));
  }
  for (int i = 1; i < intervals.size(); ++i) {
    if (intervals[i] <= 0) {
      const std::string msg("non-positive sampling interval");
      throw std::invalid_argument(msg + " -- " + std::to_string(intervals[i]));
    }
    if (intervals[i - 1] % intervals[i] != 0) {
      const std::string msg("sampling interval is no divisor of its predecessor");
      throw std::invalid_argument(msg + " -- " + std::to_string(intervals[i]));
    }
  }

  std::ofstream csv;
  if (!csvfile.empty()) {
    csv.open(csvfile + ".csv");
    if (!csv.good())
      throw std::invalid_argument("file cannot be opened -- '" + csvfile + ".csv'");
    csv << "# Input graph: " << infile << "\n";
    csv << "# OD-pairs: " << odfile << "\n";
    csv << "# Objective: " << (clp.isSet("so") ? "SO" : "UE") << "\n";
    csv << "# Function: " << clp.getValue<std::string>("f", "modified_davidson") << "\n";
    csv << "# Shortest-path algo: " << clp.getValue<std::string>("a", "dijkstra") << "\n";
    csv << "# Period of analysis: " << period << "h\n";
    csv << "# Sampling intervals: [";
    for (int i = 0, prevInterval = -1; i < intervals.size(); prevInterval = intervals[i++])
      if (intervals[i] != prevInterval)
        csv << intervals[i] << '@' << i + 1 << ':';
    csv << "1@" << intervals.size() << "]\n";
    csv << std::flush;
  }

  std::ofstream patternFile;
  if (!outfile.empty()) {
    patternFile.open(outfile + ".fp.bin", std::ios::binary);
    if (!patternFile.good())
      throw std::invalid_argument("file cannot be opened -- '" + outfile + ".fp.bin'");
    bio::write(patternFile, period);
  }

  FrankWolfeAssignmentT assign(graph, odPairs, csv, patternFile, clp.isSet("v"));

  if (csv.is_open()) {
    csv << "# Preprocessing time: " << assign.stats.totalRunningTime << "ms\n";
    csv << "iteration,sampling_interval,customization_time,query_time,line_search_time,total_time,";
    csv << "avg_change,max_change,obj_function_value,total_travel_cost,checksum\n";
    csv << std::flush;
  }

  assign.run(numIterations, intervals);
}

// Picks the shortest-path algorithm according to the command line options.
template <template <typename> class ObjFunctionT, template <typename> class TravelCostFunctionT>
void chooseShortestPathAlgo(const CommandLineParser& clp) {
  using VertexAttributes = VertexAttrs<LatLngAttribute>;
  using EdgeAttributes = EdgeAttrs<
      CapacityAttribute, EdgeIdAttribute, LengthAttribute, TravelCostAttribute,
      TravelTimeAttribute>;
  using Graph = StaticGraph<VertexAttributes, EdgeAttributes>;

  const std::string algo = clp.getValue<std::string>("a", "dijkstra");
  if (algo == "dijkstra") {
    using Assignment = FrankWolfeAssignment<
        ObjFunctionT, TravelCostFunctionT, trafficassignment::DijkstraAdapter, Graph>;
    assignTraffic<Assignment>(clp);
  } else if (algo == "bidijkstra") {
    using Assignment = FrankWolfeAssignment<
        ObjFunctionT, TravelCostFunctionT, trafficassignment::BiDijkstraAdapter, Graph>;
    assignTraffic<Assignment>(clp);
  } else if (algo == "ch") {
    using Assignment = FrankWolfeAssignment<
        ObjFunctionT, TravelCostFunctionT, trafficassignment::CHAdapter, Graph>;
    assignTraffic<Assignment>(clp);
  } else if (algo == "cch") {
    using Assignment = FrankWolfeAssignment<
        ObjFunctionT, TravelCostFunctionT, trafficassignment::CCHAdapter, Graph>;
    assignTraffic<Assignment>(clp);
  } else {
    throw std::invalid_argument("unrecognized shortest-path algorithm -- '" + algo + "'");
  }
}

// Picks the travel cost function according to the command line options.
template <template <typename> class ObjFunctionT>
void chooseTravelCostFunction(const CommandLineParser& clp) {
  const std::string func = clp.getValue<std::string>("f", "modified_davidson");
  if (func == "bpr")
    chooseShortestPathAlgo<ObjFunctionT, BprFunction>(clp);
  else if (func == "davidson")
    chooseShortestPathAlgo<ObjFunctionT, DavidsonFunction>(clp);
  else if (func == "modified_davidson")
    chooseShortestPathAlgo<ObjFunctionT, ModifiedDavidsonFunction>(clp);
  else if (func == "inverse")
    chooseShortestPathAlgo<ObjFunctionT, InverseFunction>(clp);
  else
    throw std::invalid_argument("unrecognized travel cost function -- '" + func + "'");
}

// Picks the objective function according to the command line options.
void chooseObjFunction(const CommandLineParser& clp) {
  if (clp.isSet("so"))
    chooseTravelCostFunction<SystemOptimum>(clp);
  else
    chooseTravelCostFunction<UserEquilibrium>(clp);
}

int main(int argc, char* argv[]) {
  try {
    CommandLineParser clp(argc, argv);
    if (clp.isSet("help")) {
      printUsage();
      return EXIT_SUCCESS;
    }
    chooseObjFunction(clp);
  } catch (std::invalid_argument& e) {
    std::cerr << argv[0] << ": " << e.what() << std::endl;
    std::cerr << "Try '" << argv[0] <<" -help' for more information." << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
