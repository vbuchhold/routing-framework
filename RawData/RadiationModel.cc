#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <boost/geometry.hpp>
#include <csv.h>

#include "DataStructures/Geometry/SummedAreaTables/OctagonalSummedAreaTable.h"
#include "DataStructures/Geometry/Area.h"
#include "DataStructures/Geometry/CoordinateTransformation.h"
#include "DataStructures/Geometry/Point.h"
#include "DataStructures/Graph/Attributes/LatLngAttribute.h"
#include "DataStructures/Graph/Attributes/SequentialVertexIdAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Utilities/Matrix.h"
#include "DataStructures/Utilities/OriginDestination.h"
#include "Tools/CommandLine/CommandLineParser.h"
#include "Tools/CommandLine/ProgressBar.h"
#include "Tools/Constants.h"
#include "Tools/Math.h"
#include "Tools/Timer.h"

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

void printUsage() {
  std::cout <<
      "Usage: RadiationModel [-s <factor>] -g <file> -r <file> -grid <file> -o <file>\n"
      "Generates travel demand data for the specified road network according to the\n"
      "radiation model, using the German Census 2011 population grid.\n"
      "  -s <factor>       scale the population by <factor> (defaults to 1.0)\n"
      "  -d <meters>       the max. distance between a cell's center and mapped segment\n"
      "  -seed <seed>      start the random number generator with <seed>\n"
      "  -g <file>         the network under study in binary format\n"
      "  -r <file>         restrict origins/destinations to area from OSM POLY file\n"
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

    const auto areaFilename = clp.getValue<std::string>("r");
    const auto gridFilename = clp.getValue<std::string>("grid");
    const auto graphFilename = clp.getValue<std::string>("g");
    const auto outputFilename = clp.getValue<std::string>("o");
    const double scale = clp.getValue<double>("s", 1.0);
    const int maxDist = clp.getValue<int>("d", 5000);

    // Read the graph from file.
    std::ifstream graphFile(graphFilename, std::ios::binary);
    if (!graphFile.good())
      throw std::invalid_argument("file not found -- '" + graphFilename + "'");
    StaticGraph<VertexAttrs<LatLngAttribute, SequentialVertexIdAttribute>> graph(graphFile);
    graphFile.close();

    if (graph.numVertices() > 0 && graph.sequentialVertexId(0) == INVALID_VERTEX)
      FORALL_VERTICES(graph, v)
        graph.sequentialVertexId(v) = v;

    // Read the population grid from file.
    std::vector<std::pair<::Point, int>> gridCells;
    ::Point minCell(INFTY, INFTY);
    ::Point maxCell(-INFTY, -INFTY);
    ::Point cell;
    int numInhabitants;
    io::CSVReader<3, io::trim_chars<>, io::no_quote_escape<';'>> gridFile(gridFilename);
    gridFile.read_header(io::ignore_extra_column, "x_mp_100m", "y_mp_100m", "Einwohner");
    while (gridFile.read_row(cell.x(), cell.y(), numInhabitants))
      if (numInhabitants > 0) {
        gridCells.emplace_back(cell, numInhabitants);
        minCell.min(cell);
        maxCell.max(cell);
      }

    // Allocate and fill the population matrix.
    const auto gridBounds = maxCell - minCell;
    const int numRows = gridBounds.y() / 100 + 1;
    const int numCols = gridBounds.x() / 100 + 1;
    Matrix<int> populationGrid(numRows, numCols, 0);
    for (const auto& cell : gridCells) {
      const auto pos = cell.first - minCell;
      populationGrid(numRows - pos.y() / 100 - 1, pos.x() / 100) = cell.second;
    }

    // Build an R-tree storing the road segments.
    using Point = bg::model::point<double, 2, bg::cs::spherical_equatorial<bg::degree>>;
    using Segment = bg::model::segment<Point>;
    using RoadSegment = std::tuple<Segment, int, int>;
    std::vector<RoadSegment> roadSegments;
    FORALL_VALID_EDGES(graph, u, e) {
      const int v = graph.edgeHead(e);
      const Point tail(graph.latLng(u).lngInDeg(), graph.latLng(u).latInDeg());
      const Point head(graph.latLng(v).lngInDeg(), graph.latLng(v).latInDeg());
      roadSegments.emplace_back(Segment(tail, head), u, v);
    }
    bgi::rtree<RoadSegment, bgi::quadratic<16>> rTree(roadSegments);

    // Map each cell C inside the area under study to the road segment nearest to the center of C.
    const auto primaryCrs = CoordinateTransformation::WGS_84;
    const auto secondaryCrs = CoordinateTransformation::ETRS89_LAEA_EUROPE;
    CoordinateTransformation trans(primaryCrs, secondaryCrs);
    Area area;
    area.importFromOsmPolyFile(areaFilename);
    const auto box = area.boundingBox();
    LatLng nw(box.northEast().y(), box.southWest().x());
    LatLng se(box.southWest().y(), box.northEast().x());
    double easting, northing;
    trans.forward(toRadians(nw.lngInDeg()), toRadians(nw.latInDeg()), easting, northing);
    const auto min = ::Point(std::round(easting), std::round(northing)) - minCell;
    trans.forward(toRadians(se.lngInDeg()), toRadians(se.latInDeg()), easting, northing);
    const auto max = ::Point(std::round(easting), std::round(northing)) - minCell;
    const ::Point minCoveredCell((min.x() + 50) / 100, numRows - (min.y() + 50) / 100 - 1);
    const ::Point maxCoveredCell((max.x() + 50) / 100, numRows - (max.y() + 50) / 100 - 1);
    Matrix<int> representative(numRows, numCols, INVALID_VERTEX);
    int population = 0;
    int numUnmappedCells = 0;
    for (int x = minCoveredCell.x(); x <= maxCoveredCell.x(); ++x)
      for (int y = minCoveredCell.y(); y <= maxCoveredCell.y(); ++y) {
        const auto c = ::Point(x * 100, (numRows - y - 1) * 100) + minCell;
        double lng, lat;
        trans.reverse(c.x(), c.y(), lng, lat);
        const LatLng center(toDegrees(lat), toDegrees(lng));
        if (area.contains({center.longitude(), center.latitude()})) {
          population += populationGrid(y, x);
          const Point queryPoint(center.lngInDeg(), center.latInDeg());
          RoadSegment nearestRoadSegment;
          rTree.query(bgi::nearest(queryPoint, 1), &nearestRoadSegment);
          Segment segment;
          int u, v;
          std::tie(segment, u, v) = nearestRoadSegment;
          if (bg::distance(queryPoint, segment) * EARTH_RADIUS <= maxDist) {
            const auto distToU = bg::comparable_distance(queryPoint, segment.first);
            const auto distToV = bg::comparable_distance(queryPoint, segment.second);
            representative(y, x) = distToU < distToV && graph.containsEdge(v, u) ? u : v;
          } else {
            ++numUnmappedCells;
          }
        }
      }

    // Open the output file.
    std::ofstream out(outputFilename + ".csv");
    if (!out.good())
      throw std::invalid_argument("file cannot be opened -- '" + outputFilename + ".csv'");
    out << "# Input graph: " << graphFilename << '\n';
    out << "# Methodology: radiation model (" << scale << ")\n";
    out << "origin,destination\n";

    // Generate trips between the grid cells.
    Timer timer;
    std::cout << "Generating trips: ";
    const auto coveredGridBounds = maxCoveredCell - minCoveredCell;
    ProgressBar bar((coveredGridBounds.x() + 1) * (coveredGridBounds.y() + 1));
    OctagonalSummedAreaTable sat(populationGrid);
    std::vector<OriginDestination> result;
    std::default_random_engine rand(clp.getValue<int>("seed", 19900325));
    for (int srcX = minCoveredCell.x(); srcX <= maxCoveredCell.x(); ++srcX)
      for (int srcY = minCoveredCell.y(); srcY <= maxCoveredCell.y(); ++srcY) {
        if (representative(srcY, srcX) != INVALID_VERTEX && populationGrid(srcY, srcX) > 0)
          for (int dstX = minCoveredCell.x(); dstX <= maxCoveredCell.x(); ++dstX)
            for (int dstY = minCoveredCell.y(); dstY <= maxCoveredCell.y(); ++dstY)
              if (representative(dstY, dstX) != INVALID_VERTEX && populationGrid(dstY, dstX) > 0) {
                if (srcX == dstX && srcY == dstY)
                  continue;
                const ::Point src(srcX, srcY);
                const ::Point dst(dstX, dstY);
                const double radius = src.getEuclideanDistanceTo(dst);
                int64_t srcPop = populationGrid(srcY, srcX);
                int64_t dstPop = populationGrid(dstY, dstX);
                int64_t surroundingPop = sat.sumOverOctagon(src, radius);

                // Does the surrounding population include the destination population?
                const int height = std::round(0.92387953251128675613 * radius); // r * cos(pi / 8)
                const int side = std::round(0.38268343236508977173 * radius);   // r * sin(pi / 8)
                if (src.getManhattanDistanceTo(dst) <= height + side &&
                    src.getChebyshevDistanceTo(dst) <= height)
                  surroundingPop -= dstPop;
                assert(surroundingPop > 0);

                // Pick the number of trips between the source and destination cell.
                double p = 1.0 * srcPop * dstPop / (surroundingPop * (dstPop + surroundingPop));
                assert(p > 0); assert(p <= 1);
                int numTrips = std::binomial_distribution<>(std::round(scale * srcPop), p)(rand);

                // Generate the chosen number of trips.
                for (int i = 0; i < numTrips; ++i)
                  result.emplace_back(representative(srcY, srcX), representative(dstY, dstX));
              }
        ++bar;
      }
    const int elapsed = timer.elapsed();
    std::cout << "done (" << elapsed << "ms)." << std::endl;

    std::cout << "Writing OD-pairs to file..." << std::flush;
    for (const auto& record : result) {
      out << graph.sequentialVertexId(record.origin) << ',';
      out << graph.sequentialVertexId(record.destination) << '\n';
    }
    std::cout << " done." << std::endl;
  } catch (std::exception& e) {
    std::cerr << argv[0] << ": " << e.what() << std::endl;
    std::cerr << "Try '" << argv[0] <<" -help' for more information." << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
