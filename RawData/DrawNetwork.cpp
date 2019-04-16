#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/dynamic_bitset.hpp>
#include <csv.h>

#include "Algorithms/GraphTraversal/StronglyConnectedComponents.h"
#include "DataStructures/Geometry/Area.h"
#include "DataStructures/Geometry/LatLng.h"
#include "DataStructures/Geometry/Point.h"
#include "DataStructures/Geometry/Polygon.h"
#include "DataStructures/Geometry/Rectangle.h"
#include "DataStructures/Graph/Attributes/EdgeIdAttribute.h"
#include "DataStructures/Graph/Attributes/LatLngAttribute.h"
#include "DataStructures/Graph/Attributes/NumLanesAttribute.h"
#include "DataStructures/Graph/Attributes/RoadGeometryAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Utilities/OriginDestination.h"
#include "Tools/CommandLine/CommandLineParser.h"
#include "Visualization/Graphics/PdfGraphic.h"
#include "Visualization/Graphics/PngGraphic.h"
#include "Visualization/Graphics/SvgGraphic.h"
#include "Visualization/Color.h"
#include "Visualization/Graphic.h"
#include "Visualization/PrimitiveDrawer.h"

inline void printUsage() {
  std::cout <<
      "Usage: DrawNetwork [-c <file>] -o <file> -g <file>\n"
      "       DrawNetwork [-c <file>] -o <file> -g <file> -b <file>\n"
      "       DrawNetwork [-c <file>] -o <file> -g <file> -b <file> -d <file>\n"
      "       DrawNetwork [-c <file>] -o <file> -g <file> -f <file>\n"
      "Visualizes networks, flow patterns throughout networks and travel demand data.\n"
      "  -stuttgart        remove outliers in the network of Stuttgart\n"
      "  -n                hide the network\n"
      "  -i                draw all intermediate flow patterns\n"
      "  -w <cm>           width in centimeters of the graphic (defaults to 14.0)\n"
      "  -h <cm>           height in centimeters of the graphic (defaults to 14.0)\n"
      "  -fmt <fmt>        file format of the graphic\n"
      "                      possible values: PDF PNG (default) SVG\n"
      "  -c <file>         clip the graphic to the specified OSM POLY file\n"
      "  -g <file>         draw the network in <file>\n"
      "  -b <file>         draw the boundaries in the specified OSM POLY file\n"
      "  -d <file>         draw the travel demand in <file>\n"
      "  -f <file>         draw the flow patterns in <file>\n"
      "  -o <file>         place output in <file>\n"
      "  -help             display this help and exit\n";
}

// A graph type that encompasses all attributes required for drawing.
using VertexAttributes = VertexAttrs<LatLngAttribute>;
using EdgeAttributes = EdgeAttrs<EdgeIdAttribute, NumLanesAttribute, RoadGeometryAttribute>;
using GraphT = StaticGraph<VertexAttributes, EdgeAttributes>;

// Draws the specified edge using the given primitive drawer.
inline void drawEdge(PrimitiveDrawer& pd, double width, const GraphT& graph, int u, int e) {
  pd.setLineWidth(graph.numLanes(e) * width);
  const auto v = graph.edgeHead(e);
  const auto geometry = graph.roadGeometry(e);
  if (geometry.empty()) {
    pd.drawLine(graph.latLng(u).webMercatorProjection(), graph.latLng(v).webMercatorProjection());
  } else {
    pd.drawLine(graph.latLng(u).webMercatorProjection(), geometry.front().webMercatorProjection());
    for (auto i = 1; i < geometry.size(); ++i)
      pd.drawLine(geometry[i - 1].webMercatorProjection(), geometry[i].webMercatorProjection());
    pd.drawLine(geometry.back().webMercatorProjection(), graph.latLng(v).webMercatorProjection());
  }
}

int main(int argc, char* argv[]) {
  try {
    CommandLineParser clp(argc, argv);
    if (clp.isSet("help")) {
      printUsage();
      return EXIT_SUCCESS;
    }

    // Parse the command-line options.
    const auto isStuttgartGraph = clp.isSet("stuttgart");
    const auto hideNetwork = clp.isSet("n");
    const auto drawIntermediates = clp.isSet("i");
    const auto width = clp.getValue<double>("w", 14.0);
    const auto height = clp.getValue<double>("h", 14.0);
    const auto format = clp.getValue<std::string>("fmt", "PNG");
    const auto viewportFileName = clp.getValue<std::string>("c");
    const auto graphFileName = clp.getValue<std::string>("g");
    const auto boundFileName = clp.getValue<std::string>("b");
    const auto demandFileName = clp.getValue<std::string>("d");
    const auto flowFileName = clp.getValue<std::string>("f");
    const auto outputFileName = clp.getValue<std::string>("o");

    // Read the network from file.
    std::cout << "Reading network from file..." << std::flush;
    std::ifstream graphFile(graphFileName, std::ios::binary);
    if (!graphFile.good())
      throw std::invalid_argument("file not found -- '" + graphFileName + "'");
    GraphT graph(graphFile);
    graphFile.close();
    auto id = 0;
    FORALL_EDGES(graph, e)
      graph.edgeId(e) = id++;
    const auto numEdges = graph.numEdges();
    std::vector<Point> origCoordinates(graph.numVertices());
    FORALL_VERTICES(graph, u)
      origCoordinates[u] = graph.latLng(u).webMercatorProjection();
    std::cout << " done.\n";

    if (isStuttgartGraph) {
      // Cut off the highways to Basle, Frankfurt, Zurich, Nuremberg, and Munich.
      if (graph.numVertices() != 134663 || numEdges != 307759)
        throw std::invalid_argument("unrecognized Stuttgart network");
      boost::dynamic_bitset<> bitmask(graph.numVertices());
      bitmask.set();
      bitmask[121490] = false;
      bitmask[121491] = false;
      bitmask[121492] = false;
      bitmask[121494] = false;
      bitmask[121510] = false;
      graph.extractVertexInducedSubgraph(bitmask);
      StronglyConnectedComponents scc;
      scc.run(graph);
      graph.extractVertexInducedSubgraph(scc.getLargestSccAsBitmask());
    }

    // Compute the bounding box to which the graphic is clipped.
    Rectangle boundingBox;
    if (viewportFileName.empty()) {
      FORALL_VERTICES(graph, u)
        boundingBox.extend(graph.latLng(u).webMercatorProjection());
    } else {
      Area viewport;
      viewport.importFromOsmPolyFile(viewportFileName);
      const auto box = viewport.boundingBox();
      boundingBox.extend(LatLng(box.southWest().y(), box.southWest().x()).webMercatorProjection());
      boundingBox.extend(LatLng(box.northEast().y(), box.northEast().x()).webMercatorProjection());
    }

    // Construct a graphic of the required type.
    std::unique_ptr<Graphic> graphic;
    if (format == "PDF")
      graphic.reset(new PdfGraphic(outputFileName, width, height, boundingBox));
    else if (format == "PNG")
      graphic.reset(new PngGraphic(outputFileName, width, height, boundingBox));
    else if (format == "SVG")
      graphic.reset(new SvgGraphic(outputFileName, width, height, boundingBox));
    else
      throw std::invalid_argument("unrecognized file format -- '" + format + "'");

    PrimitiveDrawer pd(graphic.get());
    if (flowFileName.empty()) {
      if (!hideNetwork) {
        // Draw the network.
        std::cout << "Drawing network..." << std::flush;
        if (!boundFileName.empty() || !demandFileName.empty())
          pd.setColor(KIT_BLACK_15);
        FORALL_VALID_EDGES(graph, u, e)
          drawEdge(pd, LineWidth::VERY_THIN, graph, u, e);
        pd.setLineWidth(LineWidth::THIN);
        std::cout << " done.\n";
      }

      if (!boundFileName.empty()) {
        // Draw the boundaries.
        std::cout << "Drawing boundaries..." << std::flush;
        Area bound;
        bound.importFromOsmPolyFile(boundFileName);
        pd.setColor(KIT_BLACK);
        for (const auto& face : bound) {
          Polygon polygon;
          for (const auto& vertex : face)
            polygon.add(LatLng(vertex.y(), vertex.x()).webMercatorProjection());
          pd.drawPolygon(polygon);
        }
        std::cout << " done.\n";
      }

      if (!demandFileName.empty()) {
        // Draw the travel demand data, each OD pair as a straight line.
        std::cout << "Drawing travel demand..." << std::flush;
        const auto odPairs = importODPairsFrom(demandFileName);
        const auto step = std::round(odPairs.size() / 131072.0);
        pd.setColor({KIT_GREEN.red(), KIT_GREEN.green(), KIT_GREEN.blue(), 3});
        for (auto i = 0; i < odPairs.size(); i += step)
          pd.drawLine(origCoordinates[odPairs[i].origin], origCoordinates[odPairs[i].destination]);
        std::cout << " done.\n";
      }
    } else {
      // Read the flow patterns from file.
      std::cout << "Reading flow patterns from file..." << std::flush;
      std::vector<double> edgeFlows;
      int iteration;
      double sat;
      int prevIteration = 1;
      io::CSVReader<
          2, io::trim_chars<>, io::no_quote_escape<','>, io::throw_on_overflow,
          io::single_line_comment<'#'>> flowFile(flowFileName);
      flowFile.read_header(io::ignore_extra_column, "iteration", "sat");
      while (flowFile.read_row(iteration, sat)) {
        if (iteration <= 0 || sat < 0)
          throw std::invalid_argument("flow file corrupt");
        if (iteration != prevIteration) {
          if (iteration < prevIteration || edgeFlows.size() % numEdges != 0)
            throw std::invalid_argument("flow file corrupt");
          prevIteration = iteration;
        }
        edgeFlows.push_back(sat);
      }
      if (edgeFlows.size() % numEdges != 0)
        throw std::invalid_argument("flow file corrupt");
      std::cout << " done.\n";

      // Draw the flow patterns, each on a distinct graphic.
      std::array<std::vector<std::pair<int, int>>, REDS_9CLASS.size() - 1> congestionLevels;
      for (auto firstEdge = 0; firstEdge < edgeFlows.size(); firstEdge += numEdges) {
        if (!(drawIntermediates || firstEdge == 0 || firstEdge == edgeFlows.size() - numEdges))
          continue;
        std::cout << "Drawing flow pattern..." << std::flush;
        if (firstEdge != 0)
          graphic->newPage();
        for (auto& level : congestionLevels)
          level.clear();
        FORALL_VALID_EDGES(graph, u, e) {
          size_t level = 100 / 20 * edgeFlows[firstEdge + graph.edgeId(e)];
          congestionLevels[std::min(level, REDS_9CLASS.size() - 2)].emplace_back(u, e);
        }
        for (int j = 0; j < congestionLevels.size(); ++j) {
          pd.setColor(REDS_9CLASS[j + 1]);
          for (const auto& edge : congestionLevels[j])
            drawEdge(pd, LineWidth::THIN, graph, edge.first, edge.second);
        }
        std::cout << " done.\n";
      }
    }
  } catch (std::exception& e) {
    std::cerr << argv[0] << ": " << e.what() << '\n';
    std::cerr << "Try '" << argv[0] <<" -help' for more information.\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
