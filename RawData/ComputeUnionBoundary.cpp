#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_set_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/enum.h>
#include <CGAL/number_utils.h>

#include "DataStructures/Geometry/Area.h"
#include "DataStructures/Geometry/CoordinateTransformation.h"
#include "DataStructures/Geometry/LatLng.h"
#include "DataStructures/Geometry/Point.h"
#include "RawData/Visum/ZonePolygons.h"
#include "Tools/CommandLine/CommandLineParser.h"
#include "Tools/Math.h"
#include "Tools/StringHelpers.h"

inline void printUsage() {
  std::cout <<
      "Usage: ComputeUnionBoundary [-col <name> -val <values>] -i <file> -o <file>\n"
      "Filters traffic zones by column, computes the union of the polygons of the\n"
      "filtered zones, and outputs the union boundary as an OSM POLY file.\n"
      "  -p <factor>       multiply coordinates in files by <factor> (defaults to 1)\n"
      "  -crs <epsg>       coordinate reference system used in the Visum files\n"
      "  -col <name>       apply filter to column <name>\n"
      "  -val <values>     space-separated list of values to be filtered for\n"
      "  -i <file>         directory that contains the Visum files\n"
      "  -o <file>         place output in <file>\n";
}

// Some CGAL primitives.
using CgalKernel = CGAL::Exact_predicates_exact_constructions_kernel;
using CgalPolygon = CGAL::Polygon_2<CgalKernel>;
using CgalPolygonWithHoles = CGAL::Polygon_with_holes_2<CgalKernel>;
using CgalPolygonSet = CGAL::Polygon_set_2<CgalKernel>;

// Converts the specified area to a CGAL polygon set.
inline CgalPolygonSet areaToCgalPolygonSet(const Area& area) {
  CgalPolygonSet cgalPolygonSet;
  CgalPolygonWithHoles cgalPolygonWithHoles;
  for (const auto& face : area) {
    if (face.orientation() == 1) {
      // Positive face, insert previous polygon with holes...
      if (!cgalPolygonWithHoles.is_unbounded()) {
        cgalPolygonSet.join(cgalPolygonWithHoles);
        cgalPolygonWithHoles.clear();
      }
      // ... and construct new outer boundary.
      for (const auto& vertex : face)
        cgalPolygonWithHoles.outer_boundary().push_back({vertex.x(), vertex.y()});
    } else {
      // Negative face, construct new hole.
      cgalPolygonWithHoles.add_hole({});
      for (const auto& vertex : face)
        (--cgalPolygonWithHoles.holes_end())->push_back({vertex.x(), vertex.y()});
    }
  }
  if (!cgalPolygonWithHoles.is_unbounded())
    cgalPolygonSet.join(cgalPolygonWithHoles);
  assert(cgalPolygonSet.is_valid());
  return cgalPolygonSet;
}

// Converts the specified CGAL polygon set to an area.
inline Area cgalPolygonSetToArea(const CgalPolygonSet& cgalPolygonSet) {
  Area area;
  std::vector<CgalPolygonWithHoles> cgalPolygonsWithHoles;
  cgalPolygonSet.polygons_with_holes(std::back_inserter(cgalPolygonsWithHoles));
  for (const auto& cgalPolygonWithHoles : cgalPolygonsWithHoles) {
    const auto& cgalOuterBoundary = cgalPolygonWithHoles.outer_boundary();
    Polygon outerBoundary;
    for (auto i = 0; i < cgalOuterBoundary.size(); ++i) {
      const auto x = std::round(CGAL::to_double(cgalOuterBoundary[i].x()));
      const auto y = std::round(CGAL::to_double(cgalOuterBoundary[i].y()));
      const Point vertex(x, y);
      if (outerBoundary.empty() || vertex != outerBoundary.back())
        outerBoundary.add(vertex);
    }
    if (outerBoundary.front() == outerBoundary.back())
      outerBoundary.removeBack();
    area.combine(outerBoundary);
    for (auto iter = cgalPolygonWithHoles.holes_begin();
         iter != cgalPolygonWithHoles.holes_end();
         ++iter) {
      const auto& cgalHole = *iter;
      Polygon hole;
      for (auto i = 0; i < cgalHole.size(); ++i) {
        const auto x = std::round(CGAL::to_double(cgalHole[i].x()));
        const auto y = std::round(CGAL::to_double(cgalHole[i].y()));
        const Point vertex(x, y);
        if (hole.empty() || vertex != hole.back())
          hole.add(Point(x, y));
      }
      if (hole.front() == hole.back())
        hole.removeBack();
      area.combine(hole);
    }
  }
  return area;
}

int main(int argc, char* argv[]) {
  try {
    CommandLineParser clp(argc, argv);
    if (clp.isSet("help")) {
      printUsage();
      return EXIT_SUCCESS;
    }

    // Parse the command-line options.
    const auto precision = clp.getValue<double>("p", 1);
    const auto crs = clp.getValue<int>("crs");
    const auto col = clp.getValue<std::string>("col", "CODE");
    const auto permittedVals = clp.getValues<std::string>("val");
    const auto infile = clp.getValue<std::string>("i");
    auto outfile = clp.getValue<std::string>("o");
    if (!endsWith(outfile, ".poly"))
      outfile += ".poly";

    // Filter the traffic zones and compute the union of the polygons of the filtered zones.
    const auto zonePolygons = visum::readZonePolygonsFrom(infile, precision, col, permittedVals);
    std::cout << "Computing union of zone polygons..." << std::flush;
    CgalPolygonSet cgalUnion;
    std::vector<CgalPolygonWithHoles> cgalPolygonsWithHoles;
    for (const auto& zone : zonePolygons) {
      cgalUnion.join(areaToCgalPolygonSet(zone.second));
      // Remove all holes from the union.
      cgalUnion.polygons_with_holes(std::back_inserter(cgalPolygonsWithHoles));
      for (auto& cgalPolygonWithHoles : cgalPolygonsWithHoles)
        for (auto iter = cgalPolygonWithHoles.holes_begin();
             iter != cgalPolygonWithHoles.holes_end();
             ++iter) {
          if (iter->is_simple()) {
            iter->reverse_orientation();
            cgalUnion.join(*iter);
          } else {
            CgalPolygon::Vertex_const_iterator i, j;
            bool cycleFound = false;
            for (i = iter->vertices_begin(); i != iter->vertices_end(); ++i) {
              for (j = i + 1; j != iter->vertices_end(); ++j)
                if (*i == *j) {
                  cycleFound = true;
                  break;
                }
              if (cycleFound)
                break;
            }
            CgalPolygon cycle(i, j);
            iter->erase(i, j);
            assert(iter->orientation() != cycle.orientation());
            if (cycle.orientation() == CGAL::CLOCKWISE) {
              cycle.reverse_orientation();
              cgalUnion.join(cycle);
              std::vector<CgalPolygonWithHoles> diff;
              difference(*iter, cycle, std::back_inserter(diff));
              assert(diff.empty());
            }
            if (iter->orientation() == CGAL::CLOCKWISE) {
              iter->reverse_orientation();
              cgalUnion.join(*iter);
              std::vector<CgalPolygonWithHoles> diff;
              difference(cycle, *iter, std::back_inserter(diff));
              assert(diff.empty());
            }
          }
        }
      cgalPolygonsWithHoles.clear();
    }
    std::cout << " done.\n";

    // Extract the largest face.
    std::cout << "Extracting largest face..." << std::flush;
    cgalUnion.polygons_with_holes(std::back_inserter(cgalPolygonsWithHoles));
    auto cgalMaxArea = cgalPolygonsWithHoles[0].outer_boundary().area();
    const auto* cgalMaxAreaPolygon = &cgalPolygonsWithHoles[0];
    for (auto i = 1; i < cgalPolygonsWithHoles.size(); ++i) {
      if (cgalPolygonsWithHoles[i].outer_boundary().area() > cgalMaxArea) {
        cgalUnion.difference(*cgalMaxAreaPolygon);
        cgalMaxArea = cgalPolygonsWithHoles[i].outer_boundary().area();
        cgalMaxAreaPolygon = &cgalPolygonsWithHoles[i];
      } else {
        cgalUnion.difference(cgalPolygonsWithHoles[i]);
      }
    }
    std::cout << " done.\n";

    // Transform to WGS84 coordinates and output the union boundary as an OSM POLY file.
    std::cout << "Transforming to WGS84 coordinates..." << std::flush;
    CoordinateTransformation trans(crs, CoordinateTransformation::WGS_84);
    const auto union_ = cgalPolygonSetToArea(cgalUnion);
    assert(union_.begin() + 1 == union_.end());
    Polygon latLngFace;
    for (const auto& vertex : *union_.begin()) {
      double lng, lat;
      trans.forward(vertex.x() / precision, vertex.y() / precision, lng, lat);
      LatLng latLng(toDegrees(lat), toDegrees(lng));
      Point latLngPoint(latLng.longitude(), latLng.latitude());
      if (latLngFace.size() >= 2 && *(latLngFace.end() - 2) == latLngPoint)
        latLngFace.removeBack();
      else if (latLngFace.size() < 1 || *(latLngFace.end() - 1) != latLngPoint)
        latLngFace.add(latLngPoint);
    }
    if (latLngFace.size() >= 2 && *(latLngFace.end() - 2) == *(latLngFace.begin()))
      latLngFace.removeRange(latLngFace.end() - 2, latLngFace.end());
    else if (latLngFace.size() >= 2 && *(latLngFace.end() - 1) == *(latLngFace.begin() + 1))
      latLngFace.removeRange(latLngFace.begin(), latLngFace.begin() + 2);
    else if (latLngFace.size() >= 1 && *(latLngFace.end() - 1) == *(latLngFace.begin()))
      latLngFace.removeRange(latLngFace.end() - 1, latLngFace.end());
    Area latLngArea;
    latLngArea.combine(latLngFace);
    std::cout << " done.\n";
    std::cout << "Writing OSM POLY file..." << std::flush;
    latLngArea.exportToOsmPolyFile(outfile);
    std::cout << " done.\n";
  } catch (std::exception& e) {
    std::cerr << argv[0] << ": " << e.what() << '\n';
    std::cerr << "Try '" << argv[0] <<" -help' for more information.\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
