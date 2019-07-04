#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <stack>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <routingkit/customizable_contraction_hierarchy.h>
#include <routingkit/nested_dissection.h>

#include "Algorithms/TrafficAssignment/Adapters/BiDijkstraAdapter.h"
#include "Algorithms/TrafficAssignment/Adapters/CCHAdapter.h"
#include "Algorithms/TrafficAssignment/Adapters/CHAdapter.h"
#include "Algorithms/TrafficAssignment/Adapters/DijkstraAdapter.h"
#include "Algorithms/TrafficAssignment/ObjectiveFunctions/SystemOptimum.h"
#include "Algorithms/TrafficAssignment/ObjectiveFunctions/UserEquilibrium.h"
#include "Algorithms/TrafficAssignment/TraversalCostFunctions/BprFunction.h"
#include "Algorithms/TrafficAssignment/TraversalCostFunctions/DavidsonFunction.h"
#include "Algorithms/TrafficAssignment/TraversalCostFunctions/InverseFunction.h"
#include "Algorithms/TrafficAssignment/TraversalCostFunctions/ModifiedDavidsonFunction.h"
#include "Algorithms/TrafficAssignment/FrankWolfeAssignment.h"
#include "DataStructures/Containers/BitVector.h"
#include "DataStructures/Graph/Attributes/CapacityAttribute.h"
#include "DataStructures/Graph/Attributes/EdgeIdAttribute.h"
#include "DataStructures/Graph/Attributes/LatLngAttribute.h"
#include "DataStructures/Graph/Attributes/LengthAttribute.h"
#include "DataStructures/Graph/Attributes/SequentialVertexIdAttribute.h"
#include "DataStructures/Graph/Attributes/TravelTimeAttribute.h"
#include "DataStructures/Graph/Attributes/TraversalCostAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Utilities/OriginDestination.h"
#include "Tools/CommandLine/CommandLineParser.h"
#include "Tools/Constants.h"
#include "Tools/StringHelpers.h"

inline void printUsage() {
  std::cout <<
      "Usage: AssignTraffic [-so] [-l] [-f <func>] [-a <algo>] -g <file> -d <file>\n"
      "Assigns OD pairs onto a network using the (conjugate) Frank-Wolfe algorithm. It\n"
      "supports different objectives, traversal cost functions and shortest-path algos.\n"
      "  -so               find the system optimum (default: user equilibrium)\n"
      "  -l                use physical lengths as metric (default: travel time)\n"
      "  -i                output all intermediate flow patterns and OD distances\n"
      "  -v                display informative messages\n"
      "  -p <hrs>          period of analysis in hours (default: 1)\n"
      "  -n <num>          number of iterations (0 means to use the stopping criterion)\n"
      "  -f <func>         traversal cost function\n"
      "                      possible values: BPR (default) Davidson M-Davidson inverse\n"
      "  -a <algo>         shortest-path algorithm\n"
      "                      possible values: Dijkstra Bi-Dijkstra CH CCH (default)\n"
      "  -o <ord>          order in which the OD pairs are processed\n"
      "                      possible values: random input sorted (default)\n"
      "  -U <num>          maximum diameter of a cell (used for ordering OD pairs)\n"
      "  -g <file>         network in binary format\n"
      "  -d <file>         OD pairs to be assigned onto the network\n"
      "  -flow <file>      place the flow pattern after each iteration in <file>\n"
      "  -dist <file>      place the OD distances after each iteration in <file>\n"
      "  -stat <file>      place statistics about the execution in <file>\n"
      "  -help             display this help and exit\n";
}

// An active vertex during a DFS, i.e., a vertex that has been reached but not finished.
struct ActiveVertex {
  // Constructs an active vertex.
  ActiveVertex(const int id, const int nextUnexploredEdge)
      : id(id), nextUnexploredEdge(nextUnexploredEdge) {}

  int id;                 // The ID of the active vertex.
  int nextUnexploredEdge; // The next unexplored incident edge.
};

// Assigns origin and destination zones to OD pairs based on a partition of the elimination tree.
template <typename GraphT>
inline void assignZonesToODPairs(
    const GraphT& inputGraph, std::vector<ClusteredOriginDestination>& odPairs, const int maxDiam) {
  // Convert the input graph to RoutingKit's graph representation.
  const int numVertices = inputGraph.numVertices();
  std::vector<float> lats(numVertices);
  std::vector<float> lngs(numVertices);
  std::vector<unsigned int> tails(inputGraph.numEdges());
  std::vector<unsigned int> heads(inputGraph.numEdges());
  FORALL_VERTICES(inputGraph, u) {
    lats[u] = inputGraph.latLng(u).latInDeg();
    lngs[u] = inputGraph.latLng(u).lngInDeg();
    FORALL_INCIDENT_EDGES(inputGraph, u, e) {
      tails[e] = u;
      heads[e] = inputGraph.edgeHead(e);
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

  // Build the elimination out-tree from the elimination in-tree.
  std::vector<int> firstChild(numVertices + 1);
  std::vector<int> children(numVertices - 1);
  for (auto v = 0; v < numVertices - 1; ++v)
    ++firstChild[tree[v]];
  auto firstEdge = 0; // The index of the first edge out of the current/next vertex.
  for (auto v = 0; v <= numVertices; ++v) {
    std::swap(firstEdge, firstChild[v]);
    firstEdge += firstChild[v];
  }
  for (auto v = 0; v < numVertices - 1; ++v)
    children[firstChild[tree[v]]++] = v;
  for (auto v = numVertices - 1; v > 0; --v)
    firstChild[v] = firstChild[v - 1];
  firstChild[0] = 0;

  // Decompose the elimination tree into as few cells with bounded diameter as possible.
  BitVector isRoot(numVertices);
  std::vector<int> height(numVertices); // height[v] is the height of the subtree rooted at v.
  for (auto v = 0; v < numVertices; ++v) {
    const auto first = firstChild[v];
    const auto last = firstChild[v + 1];
    std::sort(children.begin() + first, children.begin() + last, [&](const auto u, const auto v) {
      assert(u >= 0); assert(u < height.size());
      assert(v >= 0); assert(v < height.size());
      return height[u] < height[v];
    });
    for (auto i = first; i < last; ++i)
      if (height[v] + 1 + height[children[i]] <= maxDiam)
        height[v] = 1 + height[children[i]];
      else
        isRoot[children[i]] = true;
  }

  // Number the cells in the order in which they are discovered during a DFS from the root.
  int freeCellId = 1; // The next free cell ID.
  std::vector<int> cellIds(numVertices);
  std::stack<ActiveVertex, std::vector<ActiveVertex>> activeVertices;
  activeVertices.emplace(numVertices - 1, firstChild[numVertices - 1]);
  while (!activeVertices.empty()) {
    auto &v = activeVertices.top();
    const auto head = children[v.nextUnexploredEdge];
    ++v.nextUnexploredEdge;
    cellIds[order[head]] = isRoot[head] ? freeCellId++ : cellIds[order[v.id]];
    if (v.nextUnexploredEdge == firstChild[v.id + 1])
      activeVertices.pop();
    if (firstChild[head] != firstChild[head + 1])
      activeVertices.emplace(head, firstChild[head]);
  }

  // Assign origin and destination zones to OD pairs.
  for (auto& od : odPairs) {
    od.originZone = cellIds[od.origin];
    od.destinationZone = cellIds[od.destination];
  }
}

// Assigns all OD flows onto the graph.
template <typename FWAssignmentT>
inline void assignTraffic(const CommandLineParser& clp) {
  // Parse the command-line options.
  const auto findSO = clp.isSet("so");
  const auto useLengths = clp.isSet("l");
  const auto outputIntermediates = clp.isSet("i");
  const auto verbose = clp.isSet("v");
  const auto analysisPeriod = clp.getValue<double>("p", 0);
  const auto traversalCostFunction = clp.getValue<std::string>("f", "BPR");
  const auto shortestPathAlgorithm = clp.getValue<std::string>("a", "CCH");
  const auto ord = clp.getValue<std::string>("o", "sorted");
  const auto maxDiam = clp.getValue<int>("U", 32);
  const auto graphFileName = clp.getValue<std::string>("g");
  const auto demandFileName = clp.getValue<std::string>("d");
  auto numIterations = clp.getValue<int>("n", 0);
  auto flowFileName = clp.getValue<std::string>("flow");
  auto distFileName = clp.getValue<std::string>("dist");
  auto statFileName = clp.getValue<std::string>("stat");
  if (useLengths)
    numIterations = 1;
  if (!flowFileName.empty() && !endsWith(flowFileName, ".csv"))
    flowFileName += ".csv";
  if (!distFileName.empty() && !endsWith(distFileName, ".csv"))
    distFileName += ".csv";
  if (!statFileName.empty() && !endsWith(statFileName, ".csv"))
    statFileName += ".csv";

  // Read the graph from file.
  std::ifstream graphFile(graphFileName, std::ios::binary);
  if (!graphFile.good())
    throw std::invalid_argument("file not found -- '" + graphFileName + "'");
  typename FWAssignmentT::Graph graph(graphFile);
  graphFile.close();
  FORALL_EDGES(graph, e) {
    graph.capacity(e) = std::max(std::round(analysisPeriod * graph.capacity(e)), 1.0);
    graph.edgeId(e) = e;
  }
  if (useLengths)
    FORALL_EDGES(graph, e)
      graph.travelTime(e) = graph.length(e);
  if (graph.numVertices() > 0 && graph.sequentialVertexId(0) == INVALID_VERTEX)
    FORALL_VERTICES(graph, v) {
      assert(graph.sequentialVertexId(v) == INVALID_VERTEX);
      graph.sequentialVertexId(v) = v;
    }
  auto maxOrigId = -1;
  FORALL_VERTICES(graph, u)
    maxOrigId = std::max(maxOrigId, graph.sequentialVertexId(u));
  std::vector<int> origToCurrentId(maxOrigId + 1, INVALID_VERTEX);
  FORALL_VERTICES(graph, u) {
    assert(origToCurrentId[graph.sequentialVertexId(u)] == INVALID_VERTEX);
    origToCurrentId[graph.sequentialVertexId(u)] = u;
  }

  // Read the OD pairs from file and reorder them if necessary.
  auto odPairs = importClusteredODPairsFrom(demandFileName);
  for (auto& pair : odPairs) {
    assert(pair.origin >= 0); assert(pair.origin <= maxOrigId);
    assert(pair.destination >= 0); assert(pair.destination <= maxOrigId);
    pair.origin = origToCurrentId[pair.origin];
    pair.destination = origToCurrentId[pair.destination];
  }
  if (ord == "random") {
    std::shuffle(odPairs.begin(), odPairs.end(), std::minstd_rand());
  } else if (ord == "sorted") {
    assignZonesToODPairs(graph, odPairs, maxDiam);
    std::sort(odPairs.begin(), odPairs.end());
  } else if (ord != "input") {
    throw std::invalid_argument("unrecognized order -- '" + ord + "'");
  }

  std::ofstream flowFile;
  if (!flowFileName.empty()) {
    flowFile.open(flowFileName);
    if (!flowFile.good())
      throw std::invalid_argument("file cannot be opened -- '" + flowFileName + "'");
    if (!statFileName.empty())
      flowFile << "# Stat file: " << statFileName << "\n";
    flowFile << "iteration,vol,sat\n";
  }

  std::ofstream distFile;
  if (!distFileName.empty()) {
    distFile.open(distFileName);
    if (!distFile.good())
      throw std::invalid_argument("file cannot be opened -- '" + distFileName + "'");
    if (!statFileName.empty())
      distFile << "# Stat file: " << statFileName << "\n";
    distFile << "iteration,traversal_cost\n";
  }

  std::ofstream statFile;
  if (!statFileName.empty()) {
    statFile.open(statFileName);
    if (!statFile.good())
      throw std::invalid_argument("file cannot be opened -- '" + statFileName + "'");
    statFile << "# Graph: " << graphFileName << "\n";
    statFile << "# Demand: " << demandFileName << "\n";
    statFile << "# Objective function: " << (findSO ? "SO" : "UE") << "\n";
    statFile << "# Traversal cost function: " << traversalCostFunction << "\n";
    statFile << "# Shortest-path algorithm: " << shortestPathAlgorithm << "\n";
    statFile << "# Period of analysis: " << analysisPeriod << "\n";
    statFile << std::flush;
  }
  FWAssignmentT fwAssignment(graph, odPairs, verbose);
  if (statFile.is_open()) {
    statFile << "# Preprocessing time: " << fwAssignment.stats.totalRunningTime << "ms\n";
    statFile << "iteration,customization_time,query_time,line_search_time,total_time,";
    statFile << "prev_total_traversal_cost,prev_relative_gap,checksum\n";
    statFile << std::flush;
  }
  fwAssignment.run(flowFile, distFile, statFile, numIterations, outputIntermediates);
}

// Picks the shortest-path algorithm according to the command line options.
template <template <typename> class ObjFunctionT, template <typename> class TraversalCostFunctionT>
void chooseShortestPathAlgo(const CommandLineParser& clp) {
  using VertexAttributes = VertexAttrs<LatLngAttribute, SequentialVertexIdAttribute>;
  using EdgeAttributes = EdgeAttrs<
      CapacityAttribute, EdgeIdAttribute, LengthAttribute, TravelTimeAttribute,
      TraversalCostAttribute>;
  using Graph = StaticGraph<VertexAttributes, EdgeAttributes>;

  const auto algo = clp.getValue<std::string>("a", "CCH");
  if (algo == "Dijkstra") {
    using FWAssignment = FrankWolfeAssignment<
        ObjFunctionT, TraversalCostFunctionT, trafficassignment::DijkstraAdapter, Graph>;
    assignTraffic<FWAssignment>(clp);
  } else if (algo == "Bi-Dijkstra") {
    using FWAssignment = FrankWolfeAssignment<
        ObjFunctionT, TraversalCostFunctionT, trafficassignment::BiDijkstraAdapter, Graph>;
    assignTraffic<FWAssignment>(clp);
  } else if (algo == "CH") {
    using FWAssignment = FrankWolfeAssignment<
        ObjFunctionT, TraversalCostFunctionT, trafficassignment::CHAdapter, Graph>;
    assignTraffic<FWAssignment>(clp);
  } else if (algo == "CCH") {
    using FWAssignment = FrankWolfeAssignment<
        ObjFunctionT, TraversalCostFunctionT, trafficassignment::CCHAdapter, Graph>;
    assignTraffic<FWAssignment>(clp);
  } else {
    throw std::invalid_argument("unrecognized shortest-path algorithm -- '" + algo + "'");
  }
}

// Picks the traversal cost function according to the command line options.
template <template <typename> class ObjFunctionT>
void chooseTraversalCostFunction(const CommandLineParser& clp) {
  const auto func = clp.getValue<std::string>("f", "BPR");
  if (func == "BPR")
    chooseShortestPathAlgo<ObjFunctionT, BprFunction>(clp);
  else if (func == "Davidson")
    chooseShortestPathAlgo<ObjFunctionT, DavidsonFunction>(clp);
  else if (func == "M-Davidson")
    chooseShortestPathAlgo<ObjFunctionT, ModifiedDavidsonFunction>(clp);
  else if (func == "inverse")
    chooseShortestPathAlgo<ObjFunctionT, InverseFunction>(clp);
  else
    throw std::invalid_argument("unrecognized traversal cost function -- '" + func + "'");
}

// Picks the objective function according to the command line options.
void chooseObjFunction(const CommandLineParser& clp) {
  if (clp.isSet("so"))
    chooseTraversalCostFunction<SystemOptimum>(clp);
  else
    chooseTraversalCostFunction<UserEquilibrium>(clp);
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
