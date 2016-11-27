#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#include "Algorithms/GraphTraversal/StronglyConnectedComponents.h"
#include "DataStructures/Graph/Attributes/CapacityAttribute.h"
#include "DataStructures/Graph/Attributes/LatLngAttribute.h"
#include "DataStructures/Graph/Attributes/LengthAttribute.h"
#include "DataStructures/Graph/Attributes/TravelTimeAttribute.h"
#include "DataStructures/Graph/Attributes/XatfRoadCategoryAttribute.h"
#include "DataStructures/Graph/Export/DefaultExporter.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Graph/Import/XatfImporter.h"
#include "Tools/CommandLineParser.h"

#include "Algorithms/GraphTraversal/DfsNumbering.h"
#include "DataStructures/Graph/Attributes/AbstractBitAttribute.h"
#include "DataStructures/Graph/Attributes/RoutingCostAttribute.h"
#include "DataStructures/Graph/Attributes/TrafficFlowAttribute.h"
#include "Tools/ProgressBar.h"
#include "Tools/Timer.h"

// A graph data structure encompassing all vertex and edge attributes available for output.
using VertexAttributes = VertexAttrs<LatLngAttribute>;
using EdgeAttributes =
    EdgeAttrs<CapacityAttribute, LengthAttribute, TravelTimeAttribute, XatfRoadCategoryAttribute>;
using GraphT = StaticGraph<VertexAttributes, EdgeAttributes>;

void printUsage() {
  std::cout <<
      "Usage: ConvertGraph -s <fmt> -d <fmt> [-c] [-scc] -a <attrs> -i <file> -o <file>\n"
      "This program converts a graph from a source file format to a destination format,\n"
      "possibly extracting the largest strongly connected component of the input graph.\n"
      "  -s <fmt>          source file format\n"
      "                      possible values: default dimacs visum xatf\n"
      "  -d <fmt>          destination file format\n"
      "                      possible values: default dimacs\n"
      "  -c                compress the output file(s), if available\n"
      "  -scc              extract the largest strongly connected component\n"
      "  -a <attrs>        blank-separated list of vertex/edge attributes to be output\n"
      "                      possible values:\n"
      "                        capacity lat_lng length travel_time xatf_road_category\n"
      "  -i <file>         input file(s) without file extension\n"
      "  -o <file>         output file(s) without file extension\n"
      "  -help             display this help and exit\n";
}

// Prints the specified error message to standard error.
void printErrorMessage(const std::string& invokedName, const std::string& msg) {
  std::cerr << invokedName << ": " << msg << std::endl;
  std::cerr << "Try '" << invokedName <<" -help' for more information." << std::endl;
}

// Imports a graph according to the input file format specified on the command line and returns it.
GraphT importGraph(const CommandLineParser& clp) {
  const std::string fmt = clp.getValue<std::string>("s");

  // Choose the appropriate importer.
  if (fmt == "xatf")
    return GraphT(clp.getValue<std::string>("i"), XatfImporter());
  else
    throw std::invalid_argument("invalid input file format -- '" + fmt + "'");
}

// Executes a graph export using the specified exporter.
template <typename ExporterT>
void doExport(const CommandLineParser& clp, const GraphT& graph, ExporterT ex) {
  // Output only those attributes specified on the command line.
  std::vector<std::string> attrsToBeOutput = clp.getValues<std::string>("a");
  for (const auto& attr : GraphT::getAttributeNames())
    if (std::find(attrsToBeOutput.begin(), attrsToBeOutput.end(), attr) == attrsToBeOutput.end())
      ex.ignoreAttribute(attr);
  graph.exportTo(clp.getValue<std::string>("o"), ex);
}

// Exports the specified graph according to the output file format specified on the command line.
void exportGraph(const CommandLineParser& clp, const GraphT& graph) {
  const std::string fmt = clp.getValue<std::string>("d");
  const bool compress = clp.isSet("c");

  // Choose the appropriate exporter.
  if (fmt == "default")
    doExport(clp, graph, DefaultExporter(compress));
  else
    throw std::invalid_argument("invalid output file format -- '" + fmt + "'");
}

int main(int argc, char* argv[]) {
  CommandLineParser clp;
  try {
    clp.parse(argc, argv);
  } catch (std::invalid_argument& e) {
    printErrorMessage(argv[0], e.what());
    return EXIT_FAILURE;
  }

  if (clp.isSet("help")) {
    printUsage();
    return EXIT_SUCCESS;
  }

  try {
    std::cout << "Reading the input file(s)..." << std::flush;
    GraphT graph = importGraph(clp);
    std::cout << " done." << std::endl;

    if (clp.isSet("scc")) {
      std::cout << "Computing strongly connected components..." << std::flush;
      StronglyConnectedComponents scc;
      scc.run(graph);
      std::cout << " done." << std::endl;

      std::cout << "Extracting the largest SCC..." << std::flush;
      graph.extractVertexInducedSubgraph(scc.getLargestSccAsBitmask());
      std::cout << " done." << std::endl;
    }
    if (clp.isSet("o")) {
      std::cout << "Writing the output file(s)..." << std::flush;
      exportGraph(clp, graph);
      std::cout << " done." << std::endl;
    }
  } catch (std::invalid_argument& e) {
    printErrorMessage(argv[0], e.what());
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
