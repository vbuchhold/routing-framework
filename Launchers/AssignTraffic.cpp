#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

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
#include "Tools/CommandLine/CommandLineParser.h"

void printUsage() {
  std::cout <<
      "Usage: AssignTraffic [-f <func>] [-a <algo>] -i <file> -od <file> [-o <file>]\n"
      "This program assigns OD-pairs onto a network using the Frank-Wolfe method. It\n"
      "supports different objectives, travel cost functions and shortest-path algos.\n"
      "  -so               find the system optimum (default: user equilibrium)\n"
      "  -v                display informative messages\n"
      "  -p <hrs>          the period of analysis in hours (default: 1)\n"
      "  -f <func>         the travel cost function\n"
      "                      possible values:\n"
      "                        bpr davidson modified_davidson (default) inverse\n"
      "  -a <algo>         the shortest-path algorithm\n"
      "                      possible values: dijkstra (default) bidijkstra ch cch\n"
      "  -ord <order>      the order of the OD-pairs\n"
      "                      possible values: random input (default) sorted\n"
      "  -s <seed>         start the random number generator with <seed>\n"
      "  -i <file>         the input graph in binary format\n"
      "  -od <file>        the OD-pairs to be assigned\n"
      "  -o <file>         the output CSV file without file extension\n"
      "  -help             display this help and exit\n";
}

// Assigns all OD-flows onto the input graph.
template <typename FrankWolfeAssignmentT>
void assignTraffic(const CommandLineParser& clp) {
  const std::string infile = clp.getValue<std::string>("i");
  const std::string odfile = clp.getValue<std::string>("od");
  const std::string csvfile = clp.getValue<std::string>("o");
  const std::string ord = clp.getValue<std::string>("ord", "input");
  const float period = clp.getValue<float>("p", 1);

  std::ifstream in(infile, std::ios::binary);
  if (!in.good())
    throw std::invalid_argument("file not found -- '" + infile + "'");
  typename FrankWolfeAssignmentT::InputGraph graph(in);
  in.close();

  int id = 0;
  FORALL_EDGES(graph, e) {
    graph.capacity(e) = std::round(period * graph.capacity(e));
    graph.edgeId(e) = id++;
  }

  std::ifstream od(odfile);
  if (!od.good())
    throw std::invalid_argument("file not found -- '" + odfile + "'");
  std::vector<ClusteredOriginDestination> odPairs = importClusteredODPairsFrom(od);
  od.close();

  if (ord == "random") {
    std::default_random_engine rand(clp.getValue<int>("s", 19900325));
    std::shuffle(odPairs.begin(), odPairs.end(), rand);
  } else if (ord == "sorted") {
    std::sort(odPairs.begin(), odPairs.end());
  } else if (ord != "input") {
    throw std::invalid_argument("invalid order -- '" + ord + "'");
  }

  std::ofstream csv;
  if (!csvfile.empty()) {
    csv.open(csvfile + ".csv");
    if (!csv.good())
      throw std::invalid_argument("file cannot be opened -- '" + csvfile + ".csv'");
    csv << "# Input graph: " << infile << "\n";
    csv << "# OD-pairs: " << odfile << "\n";
    csv << "# Objective: " << (clp.isSet("so") ? "SO" : "UE") << "\n";
    csv << "# Function: " << clp.getValue<std::string>("f", "bpr") << "\n";
    csv << "# Shortest-path algo: " << clp.getValue<std::string>("a", "dijkstra") << "\n";
    csv << "# Period of analysis: " << period << "h\n";
    csv << std::flush;
  }

  FrankWolfeAssignmentT assign(graph, odPairs, csv, clp.isSet("v"));

  if (csv.is_open()) {
    csv << "# Preprocessing time: " << assign.stats.totalRunningTime << "ms\n";
    csv << "iteration,customization_time,query_time,line_search_time,total_time,";
    csv << "avg_change,max_change,total_travel_cost,checksum\n";
    csv << std::flush;
  }

  assign.run();
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
