#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include <boost/dynamic_bitset.hpp>

#include "Algorithms/GraphTraversal/StronglyConnectedComponents.h"
#include "DataStructures/Geometry/Rectangle.h"
#include "DataStructures/Graph/Attributes/CoordinateAttribute.h"
#include "DataStructures/Graph/Attributes/NumLanesAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "Tools/CommandLine/CommandLineParser.h"
#include "Visualization/Graphics/PdfGraphic.h"
#include "Visualization/Graphics/PngGraphic.h"
#include "Visualization/Graphics/SvgGraphic.h"
#include "Visualization/Graphic.h"
#include "Visualization/PrimitiveDrawer.h"

void printUsage() {
  std::cout <<
      "Usage: DrawStuttgart [-bb <box>] -g <file> -v <file> [-p <file>] -o <file>\n"
      "This program visualizes the metropolitan area of Stuttgart and flow patterns\n"
      "throughout its network.\n"
      "  -bb <box>         bounding box of the printed region\n"
      "                      possible values: inner outer (default)\n"
      "  -f <fmt>          file format of the graphic\n"
      "                      possible values: pdf png (default) svg\n"
      "  -w <cm>           width in cm of the graphic (defaults to 14)\n"
      "  -h <cm>           height in cm of the graphic (defaults to 14)\n"
      "  -g <file>         network of Stuttgart in binary format\n"
      "  -v <file>         Visum file that contains polylines representing edges\n"
      "  -p <file>         flow pattern file resulting from traffic assignment\n"
      "  -o <file>         place the output in <file>\n"
      "  -help             display this help and exit\n";
}

int main(int argc, char* argv[]) {
  try {
    CommandLineParser clp(argc, argv);
    if (clp.isSet("help")) {
      printUsage();
      return EXIT_SUCCESS;
    }

    // Read the network of Stuttgart from disk.
    const std::string& graphFilename = clp.getValue<std::string>("g");
    std::ifstream graphFile(graphFilename, std::ios::binary);
    if (!graphFile.good())
      throw std::invalid_argument("file not found -- '" + graphFilename + "'");
    StaticGraph<VertexAttrs<CoordinateAttribute>, EdgeAttrs<NumLanesAttribute>> graph(graphFile);
    graphFile.close();
    if (graph.numVertices() != 134663)
      throw std::invalid_argument("invalid network");

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

    // Select the bounding box of the printed region.
    Rectangle boundingBox;
    const std::string& box = clp.getValue<std::string>("bb", "outer");
    if (box == "inner")
      boundingBox = {{3502933, 5394965}, {3523284, 5414346}};
    else if (box == "outer")
      boundingBox = {{3419985, 5322402}, {3606696, 5476854}};
    else
      throw std::invalid_argument("invalid bounding box -- '" + box + "'");

    // Create a graphic of the specified type.
    std::unique_ptr<Graphic> graphic;
    const std::string& fmt = clp.getValue<std::string>("f", "png");
    const std::string& outfile = clp.getValue<std::string>("o");
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
    if (!clp.isSet("p"))
      // Draw only the network, without any traffic flows.
      FORALL_VALID_EDGES(graph, u, e) {
        const double width = graph.numLanes(e) * LineWidth::VERY_THIN;
        pd.drawLine(graph.coordinate(u), graph.coordinate(graph.edgeHead(e)), width);
      }
    else
      throw std::invalid_argument("unsupported operation");
  } catch (std::exception& e) {
    std::cerr << argv[0] << ": " << e.what() << std::endl;
    std::cerr << "Try '" << argv[0] <<" -help' for more information." << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
