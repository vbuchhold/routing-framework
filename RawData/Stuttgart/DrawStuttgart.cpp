#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/dynamic_bitset.hpp>
#include <csv.h>

#include "Algorithms/GraphTraversal/StronglyConnectedComponents.h"
#include "DataStructures/Geometry/Area.h"
#include "DataStructures/Geometry/CoordinateConversion.h"
#include "DataStructures/Geometry/LatLng.h"
#include "DataStructures/Geometry/Point.h"
#include "DataStructures/Geometry/Polygon.h"
#include "DataStructures/Geometry/Rectangle.h"
#include "DataStructures/Graph/Attributes/CapacityAttribute.h"
#include "DataStructures/Graph/Attributes/CoordinateAttribute.h"
#include "DataStructures/Graph/Attributes/EdgeIdAttribute.h"
#include "DataStructures/Graph/Attributes/NumLanesAttribute.h"
#include "DataStructures/Graph/Attributes/VertexIdAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Utilities/OriginDestination.h"
#include "Tools/CommandLine/CommandLineParser.h"
#include "Tools/Constants.h"
#include "Visualization/Graphics/PdfGraphic.h"
#include "Visualization/Graphics/PngGraphic.h"
#include "Visualization/Graphics/SvgGraphic.h"
#include "Visualization/Color.h"
#include "Visualization/Graphic.h"
#include "Visualization/PrimitiveDrawer.h"

void printUsage() {
  std::cout <<
      "Usage: DrawStuttgart -g <file> -v <file> -p <file> -b <file> -d <file> -o <file>\n"
      "Visualize the metropolitan area of Stuttgart, flow patterns throughout its\n"
      "network, and travel demand data.\n"
      "  -c <box>          clip the output to <box>\n"
      "                      possible values: tiny small normal large (default)\n"
      "  -f <fmt>          file format of the graphic\n"
      "                      possible values: pdf png (default) svg\n"
      "  -w <cm>           width in cm of the graphic (defaults to 14)\n"
      "  -h <cm>           height in cm of the graphic (defaults to 14)\n"
      "  -a <hrs>          analysis period in hours (defaults to 1.0)\n"
      "  -g <file>         network of Stuttgart in binary format\n"
      "  -v <file>         Visum file that contains polylines representing edges\n"
      "  -p <file>         flow pattern file resulting from traffic assignment\n"
      "  -b <file>         boundary in OSM POLY format\n"
      "  -d <file>         travel demand file that contains the OD-pairs\n"
      "  -o <file>         place the output in <file>\n"
      "  -help             display this help and exit\n";
}

// A graph type that encompasses all attributes required for drawing.
using VertexAttributes = VertexAttrs<CoordinateAttribute, VertexIdAttribute>;
using EdgeAttributes = EdgeAttrs<CapacityAttribute, EdgeIdAttribute, NumLanesAttribute>;
using GraphT = StaticGraph<VertexAttributes, EdgeAttributes>;

// A map type for storing detailed edge polylines.
using EdgePolylineMap = std::map<std::pair<int, int>, std::vector<Point>>;

// Draws the specified edge using the given primitive drawer.
void drawEdge(
    PrimitiveDrawer& pd, const double width,
    const GraphT& graph, const EdgePolylineMap& edgePolylines, const int u, const int e) {
  pd.setLineWidth(graph.numLanes(e) * width);
  const int v = graph.edgeHead(e);
  const auto entry = edgePolylines.find(std::minmax(graph.vertexId(u), graph.vertexId(v)));
  if (entry == edgePolylines.end()) {
    pd.drawLine(graph.coordinate(u), graph.coordinate(v));
  } else {
    pd.drawLine(graph.coordinate(std::min(u, v)), entry->second.front());
    for (int i = 1; i < entry->second.size(); ++i)
      pd.drawLine(entry->second[i - 1], entry->second[i]);
    pd.drawLine(entry->second.back(), graph.coordinate(std::max(u, v)));
  }
}

int main(int argc, char* argv[]) {
  try {
    CommandLineParser clp(argc, argv);
    if (clp.isSet("help")) {
      printUsage();
      return EXIT_SUCCESS;
    }

    // Select the bounding box to be used as clipping path.
    Rectangle boundingBox;
    const std::string box = clp.getValue<std::string>("c", "large");
    if (box == "tiny")
      boundingBox = {{3508021, 5399810}, {3518196, 5409501}};
    else if (box == "small")
      boundingBox = {{3502933, 5394965}, {3523284, 5414346}};
    else if (box == "normal")
      boundingBox = {{3481370, 5373316}, {3573272, 5438087}};
    else if (box == "large")
      boundingBox = {{3419985, 5322402}, {3606696, 5476854}};
    else
      throw std::invalid_argument("invalid bounding box -- '" + box + "'");

    // Read the network of Stuttgart from file.
    const std::string graphFilename = clp.getValue<std::string>("g");
    std::ifstream graphFile(graphFilename, std::ios::binary);
    if (!graphFile.good())
      throw std::invalid_argument("file not found -- '" + graphFilename + "'");
    GraphT graph(graphFile);
    graphFile.close();
    if (graph.numVertices() != 134663)
      throw std::invalid_argument("invalid network");
    int id = 0;
    FORALL_EDGES(graph, e)
      graph.edgeId(e) = id++;
    const std::vector<Point> origCoordinates(&graph.coordinate(0), &graph.coordinate(134662) + 1);

    // Cut off the highways to Basle, Frankfurt, Zurich, Nuremberg, and Munich.
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

    // Read the edge polylines from file.
    EdgePolylineMap edgePolylines;
    std::pair<int, int> edge;
    int idx;
    Point coordinate;
    std::pair<int, int> prevEdge = {INVALID_ID, INVALID_ID};
    const std::string visumFilename = clp.getValue<std::string>("v");
    io::CSVReader<5, io::trim_chars<>, io::no_quote_escape<';'>> visumFile(visumFilename);
    const auto ignore = io::ignore_extra_column;
    visumFile.read_header(ignore, "VONKNOTNR", "NACHKNOTNR", "INDEX", "XKOORD", "YKOORD");
    while (visumFile.read_row(edge.first, edge.second, idx, coordinate.getX(), coordinate.getY())) {
      if (edge.first < 0 || edge.second < 0)
        throw std::invalid_argument("Visum file corrupt");
      if (edge != prevEdge) {
        if (idx != 1)
          throw std::invalid_argument("Visum file corrupt");
        if (prevEdge.first > prevEdge.second) {
          auto& prev = edgePolylines[{prevEdge.second, prevEdge.first}];
          std::reverse(prev.begin(), prev.end());
        }
        prevEdge = edge;
      }
      auto& polyline = edgePolylines[std::minmax(edge.first, edge.second)];
      polyline.push_back(coordinate);
      if (polyline.size() != idx)
        throw std::invalid_argument("Visum file corrupt");
    }
    if (prevEdge.first > prevEdge.second) {
      auto& prev = edgePolylines[{prevEdge.second, prevEdge.first}];
      std::reverse(prev.begin(), prev.end());
    }

    // Create a graphic of the specified type.
    std::unique_ptr<Graphic> graphic;
    const std::string fmt = clp.getValue<std::string>("f", "png");
    const std::string outfile = clp.getValue<std::string>("o");
    const double width = clp.getValue<double>("w", 14);
    const double height = clp.getValue<double>("h", 14);
    if (fmt == "pdf")
      graphic.reset(new PdfGraphic(outfile, width, height, boundingBox));
    else if (fmt == "png")
      graphic.reset(new PngGraphic(outfile, width, height, boundingBox));
    else if (fmt == "svg")
      graphic.reset(new SvgGraphic(outfile, width, height, boundingBox));
    else
      throw std::invalid_argument("invalid graphics format -- '" + fmt + "'");

    PrimitiveDrawer pd(graphic.get());
    if (!clp.isSet("p")) {
      // Draw only the network, without any traffic flows.
      if (clp.isSet("b"))
        pd.setColor({217, 217, 217});
      FORALL_VALID_EDGES(graph, u, e)
        drawEdge(pd, LineWidth::VERY_THIN, graph, edgePolylines, u, e);
      pd.setLineWidth(LineWidth::THIN);

      if (clp.isSet("b")) {
        // Draw the boundary.
        pd.setColor(KIT_BLACK);
        CoordinateConversion conv(CoordinateConversion::DHDN_GAUSS_KRUGER_ZONE_3);
        Area area;
        area.importFromOsmPolyFile(clp.getValue<std::string>("b"));
        for (auto& face : area) {
          Polygon polygon;
          for (const auto& vertex : face)
            polygon.add(conv.convert(LatLng(vertex.getY(), vertex.getX())));
          assert(polygon.simple());
          pd.drawPolygon(polygon);
        }
      }

      if (clp.isSet("d")) {
        // Draw the travel demand data, each OD-pair as a straight line.
        pd.setColor({0, 150, 130, 7});
        const auto odPairs = importODPairsFrom(clp.getValue<std::string>("d"));
        for (const auto& od : odPairs)
          pd.drawLine(origCoordinates[od.origin], origCoordinates[od.destination]);
      }
    } else {
      // Read the pattern file into memory.
      std::vector<std::vector<double>> flowPatterns;
      int iteration;
      double flow;
      int prevIteration = 0;
      io::CSVReader<
          2, io::trim_chars<>, io::no_quote_escape<','>, io::throw_on_overflow,
          io::single_line_comment<'#'>> patternFile(clp.getValue<std::string>("p"));
      patternFile.read_header(io::ignore_no_column, "iteration", "edge_flow");
      while (patternFile.read_row(iteration, flow)) {
        if (iteration <= 0 || flow < 0)
          throw std::invalid_argument("pattern file corrupt");
        if (iteration != prevIteration) {
          if (!flowPatterns.empty() && flowPatterns.back().size() != 307759)
            throw std::invalid_argument("pattern file corrupt");
          flowPatterns.emplace_back();
          ++prevIteration;
        }
        flowPatterns.back().push_back(flow);
      }

      // Scale segment capacities according to the analysis period.
      const double period = clp.getValue<double>("a", 1.0);
      FORALL_EDGES(graph, e)
        graph.capacity(e) = std::max(std::round(period * graph.capacity(e)), 1.0);

      // Draw the flow pattern after each iteration on a separate graphic.
      std::array<std::vector<std::pair<int, int>>, REDS_9CLASS.size() - 1> congestionLevels;
      for (int i = 0; i < flowPatterns.size(); ++i) {
        if (i != 0)
          graphic->newPage();
        for (auto& congestionLevel : congestionLevels)
          congestionLevel.clear();
        FORALL_VALID_EDGES(graph, u, e) {
          const size_t level = 100 / 20 * flowPatterns[i][graph.edgeId(e)] / graph.capacity(e);
          congestionLevels[std::min(level, REDS_9CLASS.size() - 2)].emplace_back(u, e);
        }
        for (int j = 0; j < congestionLevels.size(); ++j) {
          pd.setColor(REDS_9CLASS[j + 1]);
          for (const auto& edge : congestionLevels[j])
            drawEdge(pd, LineWidth::THIN, graph, edgePolylines, edge.first, edge.second);
        }
      }
    }
  } catch (std::exception& e) {
    std::cerr << argv[0] << ": " << e.what() << std::endl;
    std::cerr << "Try '" << argv[0] <<" -help' for more information." << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
