#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <csv.h>
#include <routingkit/geo_position_to_node.h>

#include "DataStructures/Geometry/SummedAreaTables/OctagonalSummedAreaTable.h"
#include "DataStructures/Geometry/Area.h"
#include "DataStructures/Geometry/CoordinateConversion.h"
#include "DataStructures/Geometry/Point.h"
#include "DataStructures/Graph/Attributes/LatLngAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Utilities/Matrix.h"
#include "Tools/CommandLine/CommandLineParser.h"
#include "Tools/CommandLine/ProgressBar.h"
#include "Tools/Constants.h"
#include "Tools/Timer.h"

void printUsage() {
  std::cout <<
      "Usage: RadiationModel [-p <prob>] -g <file> -a <file> -grid <file> -o <file>\n"
      "Generate a travel demand for the specified road network according to the\n"
      "radiation model, using the German Census 2011 population grid.\n"
      "  -p <prob>         the prob that an individual starts a journey (default: 0.45)\n"
      "  -s <seed>         the seed for the random-number generator\n"
      "  -g <file>         the network under study in binary format\n"
      "  -a <file>         restrict origins/destinations to area from OSM POLY file\n"
      "  -grid <file>      the German Census 2011 population grid\n"
      "  -o <file>         place output in <file>\n"
      "  -help             display this help and exit\n";
}

int main(int argc, char* argv[]) {
  try {
    CommandLineParser clp(argc, argv);
    if (clp.isSet("help")) {
      printUsage();
      return EXIT_SUCCESS;
    }

    const auto areaFilename = clp.getValue<std::string>("a");
    const auto gridFilename = clp.getValue<std::string>("grid");
    const auto graphFilename = clp.getValue<std::string>("g");
    const auto outputFilename = clp.getValue<std::string>("o");
    const double prob = clp.getValue<double>("p", 0.45);

    // Read the graph into memory.
    std::ifstream graphFile(graphFilename, std::ios::binary);
    if (!graphFile.good())
      throw std::invalid_argument("file not found -- '" + graphFilename + "'");
    StaticGraph<VertexAttrs<LatLngAttribute>, EdgeAttrs<>> graph(graphFile);
    graphFile.close();

    // Read the population grid into memory.
    std::vector<std::pair<Point, int>> gridCells;
    Point minCell(INFTY, INFTY);
    Point maxCell(-INFTY, -INFTY);
    Point cell;
    int numInhabitants;
    io::CSVReader<3, io::trim_chars<>, io::no_quote_escape<';'>> gridFile(gridFilename);
    gridFile.read_header(io::ignore_extra_column, "x_mp_100m", "y_mp_100m", "Einwohner");
    while (gridFile.read_row(cell.getX(), cell.getY(), numInhabitants))
      if (numInhabitants > 0) {
        gridCells.emplace_back(cell, numInhabitants);
        minCell.min(cell);
        maxCell.max(cell);
      }

    // Build a vp-tree for the vertices in the graph.
    std::vector<float> lat(graph.numVertices());
    std::vector<float> lng(graph.numVertices());
    FORALL_VERTICES(graph, v) {
      lat[v] = graph.latLng(v).latInDeg();
      lng[v] = graph.latLng(v).lngInDeg();
    }
    RoutingKit::GeoPositionToNode vpTree(lat, lng);

    // Allocate and fill the population matrix.
    auto gridBounds = maxCell - minCell;
    int numRows = gridBounds.getY() / 100 + 1;
    int numCols = gridBounds.getX() / 100 + 1;
    Matrix<int> populationGrid(numRows, numCols, 0);
    for (const auto& cell : gridCells) {
      const auto pos = cell.first - minCell;
      populationGrid(numRows - pos.getY() / 100 - 1, pos.getX() / 100) = cell.second;
    }

    // Obtain for each cell C inside the area under study the vertex closest to the centroid of C.
    CoordinateConversion conv(3035);
    Area area;
    area.importFromOsmPolyFile(areaFilename);
    auto box = area.boundingBox();
    auto min = conv.convert(LatLng(box.getNorthEast().getY(), box.getSouthWest().getX())) - minCell;
    auto max = conv.convert(LatLng(box.getSouthWest().getY(), box.getNorthEast().getX())) - minCell;
    Point minCoveredCell((min.getX() + 50) / 100, numRows - (min.getY() + 50) / 100 - 1);
    Point maxCoveredCell((max.getX() + 50) / 100, numRows - (max.getY() + 50) / 100 - 1);
    Matrix<int> representative(numRows, numCols, INVALID_VERTEX);
    int population = 0;
    for (int x = minCoveredCell.getX(); x <= maxCoveredCell.getX(); ++x)
      for (int y = minCoveredCell.getY(); y <= maxCoveredCell.getY(); ++y) {
        const auto centroid = conv.convert(Point(x * 100, (numRows - y - 1) * 100) + minCell);
        if (area.contains({centroid.longitude(), centroid.latitude()})) {
          representative(y, x) = vpTree.find_nearest_neighbor_within_radius(
              centroid.latInDeg(), centroid.lngInDeg(), std::numeric_limits<float>::max()).id;
          population += populationGrid(y, x);
        }
      }
    std::cout << "Population in area: " << population << std::endl;

    // Open the output file.
    std::ofstream out(outputFilename + ".csv");
    if (!out.good())
      throw std::invalid_argument("file cannot be opened -- '" + outputFilename + ".csv'");
    out << "# Input graph: " << graphFilename << '\n';
    out << "# Methodology: radiation model (" << prob << ")\n";
    out << "origin,destination\n";

    // Generate trips between the grid cells.
    OctagonalSummedAreaTable sat(populationGrid);
    std::default_random_engine rand(clp.getValue<int>("s", 19900325));
    const auto coveredGridBounds = maxCoveredCell - minCoveredCell;
    std::cout << "Generating trips: ";
    ProgressBar bar((coveredGridBounds.getX() + 1) * (coveredGridBounds.getY() + 1));
    Timer timer;
    for (int srcX = minCoveredCell.getX(); srcX <= maxCoveredCell.getX(); ++srcX)
      for (int srcY = minCoveredCell.getY(); srcY <= maxCoveredCell.getY(); ++srcY) {
        if (representative(srcY, srcX) != INVALID_VERTEX && populationGrid(srcY, srcX) != 0)
          for (int dstX = minCoveredCell.getX(); dstX <= maxCoveredCell.getX(); ++dstX)
            for (int dstY = minCoveredCell.getY(); dstY <= maxCoveredCell.getY(); ++dstY)
              if (representative(dstY, dstX) != INVALID_VERTEX && populationGrid(dstY, dstX) != 0) {
                if (representative(srcY, srcX) == representative(dstY, dstX))
                  continue;
                Point src(srcX, srcY);
                Point dst(dstX, dstY);
                double radius = src.getEuclideanDistanceTo(dst);
                int64_t srcPop = populationGrid(srcY, srcX);
                int64_t dstPop = populationGrid(dstY, dstX);
                int64_t surroundingPop = sat.sumOverOctagon(src, radius);

                // Does the surrounding population include the destination population?
                int height = std::round(0.92387953251128675613 * radius); // r * cos(pi / 8)
                int side = std::round(0.38268343236508977173 * radius);   // r * sin(pi / 8)
                if (src.getManhattanDistanceTo(dst) <= height + side &&
                    src.getChebyshevDistanceTo(dst) <= height)
                  surroundingPop -= dstPop;
                assert(surroundingPop > 0);

                // Pick the number of trips between the source and destination cell.
                double p = prob * srcPop * dstPop / (surroundingPop * (dstPop + surroundingPop));
                assert(p > 0); assert(p <= 1);
                int numTrips = std::binomial_distribution<>(srcPop, p)(rand);

                // Generate the chosen number of trips.
                for (int i = 0; i < numTrips; ++i)
                  out << representative(srcY, srcX) << ',' << representative(dstY, dstX) << '\n';
              }
        ++bar;
      }
    const int elapsed = timer.elapsed<std::chrono::seconds>();
    std::cout << "done (" << elapsed << "s)." << std::endl;
  } catch (std::exception& e) {
    std::cerr << argv[0] << ": " << e.what() << std::endl;
    std::cerr << "Try '" << argv[0] <<" -help' for more information." << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
