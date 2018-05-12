#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "Algorithms/GraphTraversal/StronglyConnectedComponents.h"
#include "DataStructures/Geometry/Area.h"
#include "DataStructures/Geometry/Point.h"
#include "DataStructures/Graph/Attributes/CapacityAttribute.h"
#include "DataStructures/Graph/Attributes/CoordinateAttribute.h"
#include "DataStructures/Graph/Attributes/FreeFlowSpeedAttribute.h"
#include "DataStructures/Graph/Attributes/LatLngAttribute.h"
#include "DataStructures/Graph/Attributes/LengthAttribute.h"
#include "DataStructures/Graph/Attributes/NumLanesAttribute.h"
#include "DataStructures/Graph/Attributes/OsmRoadCategoryAttribute.h"
#include "DataStructures/Graph/Attributes/RoadGeometryAttribute.h"
#include "DataStructures/Graph/Attributes/SequentialVertexIdAttribute.h"
#include "DataStructures/Graph/Attributes/SpeedLimitAttribute.h"
#include "DataStructures/Graph/Attributes/TravelTimeAttribute.h"
#include "DataStructures/Graph/Attributes/VertexIdAttribute.h"
#include "DataStructures/Graph/Attributes/XatfRoadCategoryAttribute.h"
#include "DataStructures/Graph/Export/DefaultExporter.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Graph/Import/OsmImporter.h"
#include "DataStructures/Graph/Import/VisumImporter.h"
#include "DataStructures/Graph/Import/XatfImporter.h"
#include "Tools/CommandLine/CommandLineParser.h"
#include "Tools/ContainerHelpers.h"

// A graph data structure encompassing all vertex and edge attributes available for output.
using VertexAttributes = VertexAttrs<
    CoordinateAttribute, LatLngAttribute, SequentialVertexIdAttribute, VertexIdAttribute>;
using EdgeAttributes = EdgeAttrs<
    CapacityAttribute, FreeFlowSpeedAttribute, LengthAttribute,
    NumLanesAttribute, OsmRoadCategoryAttribute, RoadGeometryAttribute, SpeedLimitAttribute,
    TravelTimeAttribute, XatfRoadCategoryAttribute>;
using GraphT = StaticGraph<VertexAttributes, EdgeAttributes>;

void printUsage() {
  std::cout <<
      "Usage: ConvertGraph -s <fmt> -d <fmt> [-c] [-scc] -a <attrs> -i <file> -o <file>\n"
      "This program converts a graph from a source file format to a destination format,\n"
      "possibly extracting the largest strongly connected component of the input graph.\n"
      "  -s <fmt>          source file format\n"
      "                      possible values: binary default dimacs osm visum xatf\n"
      "  -d <fmt>          destination file format\n"
      "                      possible values: binary default dimacs\n"
      "  -c                compress the output file(s), if available\n"
      "  -p <file>         extract a region given as an OSM POLY file\n"
      "  -scc              extract the largest strongly connected component\n"
      "  -ts <sys>         the system whose network is to be imported (Visum only)\n"
      "  -cs <epsg-code>   input coordinate system (Visum only)\n"
      "  -ap <hrs>         analysis period, capacity is in vehicles/AP (Visum only)\n"
      "  -a <attrs>        blank-separated list of vertex/edge attributes to be output\n"
      "                      possible values:\n"
      "                        capacity coordinate free_flow_speed lat_lng length\n"
      "                        num_lanes osm_road_category road_geometry\n"
      "                        sequential_vertex_id speed_limit travel_time vertex_id\n"
      "                        xatf_road_category\n"
      "  -i <file>         input file(s) without file extension\n"
      "  -o <file>         output file(s) without file extension\n"
      "  -help             display this help and exit\n";
}

// Imports a graph according to the input file format specified on the command line and returns it.
GraphT importGraph(const CommandLineParser& clp) {
  const std::string fmt = clp.getValue<std::string>("s");
  const std::string infile = clp.getValue<std::string>("i");

  // Pick the appropriate import procedure.
  if (fmt == "binary") {
    std::ifstream in(infile + ".gr.bin", std::ios::binary);
    if (!in.good())
      throw std::invalid_argument("file not found -- '" + infile + ".gr.bin'");
    return GraphT(in);
  } else if (fmt == "osm") {
    return GraphT(infile, OsmImporter());
  } else if (fmt == "visum") {
    const std::string sys = clp.getValue<std::string>("ts", "P");
    const int crs = clp.getValue<int>("cs", 31467);
    const int ap = clp.getValue<int>("ap", 24);
    if (ap <= 0) {
      const auto what = "analysis period not strictly positive -- '" + std::to_string(ap) + "'";
      throw std::invalid_argument(what);
    }
    return GraphT(infile, VisumImporter(infile, sys, crs, ap));
  } else if (fmt == "xatf") {
    return GraphT(infile, XatfImporter());
  } else {
    throw std::invalid_argument("unrecognized input file format -- '" + fmt + "'");
  }
}

// Executes a graph export using the specified exporter.
template <typename ExporterT>
void doExport(const CommandLineParser& clp, const GraphT& graph, ExporterT ex) {
  // Output only those attributes specified on the command line.
  std::vector<std::string> attrsToOutput = clp.getValues<std::string>("a");
  for (const auto& attr : GraphT::getAttributeNames())
    if (!contains(attrsToOutput.begin(), attrsToOutput.end(), attr))
      ex.ignoreAttribute(attr);
  graph.exportTo(clp.getValue<std::string>("o"), ex);
}

// Exports the specified graph according to the output file format specified on the command line.
void exportGraph(const CommandLineParser& clp, const GraphT& graph) {
  const std::string fmt = clp.getValue<std::string>("d");
  const bool compress = clp.isSet("c");

  // Pick the appropriate export procedure.
  if (fmt == "binary") {
    const std::string outfile = clp.getValue<std::string>("o");
    std::ofstream out(outfile + ".gr.bin", std::ios::binary);
    if (!out.good())
      throw std::invalid_argument("file cannot be opened -- '" + outfile + ".gr.bin'");
    // Output only those attributes specified on the command line.
    std::vector<std::string> attrsToIgnore;
    std::vector<std::string> attrsToOutput = clp.getValues<std::string>("a");
    for (const auto& attr : GraphT::getAttributeNames())
      if (!contains(attrsToOutput.begin(), attrsToOutput.end(), attr))
        attrsToIgnore.push_back(attr);
    graph.writeTo(out, attrsToIgnore);
  } else if (fmt == "default") {
    doExport(clp, graph, DefaultExporter(compress));
  } else {
    throw std::invalid_argument("unrecognized output file format -- '" + fmt + "'");
  }
}

int main(int argc, char* argv[]) {
  try {
    CommandLineParser clp(argc, argv);
    if (clp.isSet("help")) {
      printUsage();
      return EXIT_SUCCESS;
    }

    std::cout << "Reading the input file(s)..." << std::flush;
    auto graph = importGraph(clp);
    std::cout << " done." << std::endl;

    if (clp.isSet("p")) {
      std::cout << "Extracting the given region..." << std::flush;
      boost::dynamic_bitset<> isVertexInsideRegion(graph.numVertices());
      Area area;
      area.importFromOsmPolyFile(clp.getValue<std::string>("p"));
      const auto box = area.boundingBox();
      FORALL_VERTICES(graph, v) {
        const Point p(graph.latLng(v).longitude(), graph.latLng(v).latitude());
        isVertexInsideRegion[v] = box.contains(p) && area.contains(p);
      }
      FORALL_VERTICES(graph, v)
        graph.sequentialVertexId(v) = v;
      graph.extractVertexInducedSubgraph(isVertexInsideRegion);
      std::cout << " done." << std::endl;
    }

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
    std::cerr << argv[0] << ": " << e.what() << std::endl;
    std::cerr << "Try '" << argv[0] <<" -help' for more information." << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
