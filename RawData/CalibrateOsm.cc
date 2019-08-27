#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <routingkit/geo_position_to_node.h>

#include "DataStructures/Geometry/LatLng.h"
#include "DataStructures/Graph/Attributes/CapacityAttribute.h"
#include "DataStructures/Graph/Attributes/FreeFlowSpeedAttribute.h"
#include "DataStructures/Graph/Attributes/LatLngAttribute.h"
#include "DataStructures/Graph/Attributes/NumLanesAttribute.h"
#include "DataStructures/Graph/Attributes/OsmRoadCategoryAttribute.h"
#include "DataStructures/Graph/Attributes/SpeedLimitAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "Tools/CommandLine/CommandLineParser.h"
#include "Tools/Constants.h"

inline void printUsage() {
  std::cout <<
      "CalibrateOsm [-d <meters>] -osm <file> -ref <file> -o <file>\n"
      "This program takes an OSM and a reference graph as input, matches the edges in\n"
      "both graphs with each other, and outputs several attributes for each edge.\n"
      "  -d <meters>       only match vertices within a maximum distance of <meters>\n"
      "  -osm <file>       the OSM graph in binary format\n"
      "  -ref <file>       the reference graph in binary format\n"
      "  -o <file>         place output in <file>\n";
}

using VPTree = RoutingKit::GeoPositionToNode; // An implementation of a vantage-point tree.

// Returns the point in the specified vp-tree that is closest to the query point p.
inline int findNearestNeighbor(const VPTree& tree, const LatLng& p, const int maxDist) {
  auto neighbors = tree.find_all_nodes_within_radius(p.latInDeg(), p.lngInDeg(), maxDist);
  std::sort(neighbors.begin(), neighbors.end(), [](const auto& u, const auto& v) {
    return u.distance < v.distance;
  });
  const int numNeighbors = neighbors.size();
  if (numNeighbors == 1 || (numNeighbors > 1 && neighbors[0].distance < neighbors[1].distance))
    return neighbors[0].id;
  else
    return INVALID_INDEX;
}

int main(int argc, char* argv[]) {
  try {
    CommandLineParser clp(argc, argv);
    if (clp.isSet("help")) {
      printUsage();
      return EXIT_SUCCESS;
    }

    const auto osmFilename = clp.getValue<std::string>("osm");
    const auto refFilename = clp.getValue<std::string>("ref");
    const auto outfile = clp.getValue<std::string>("o");
    const int maxDist = clp.getValue<int>("d", 30);

    // Read both input graphs from disk.
    using VertexAttributes = VertexAttrs<LatLngAttribute>;
    using EdgeAttributes = EdgeAttrs<
        CapacityAttribute, FreeFlowSpeedAttribute,
        NumLanesAttribute, OsmRoadCategoryAttribute, SpeedLimitAttribute>;
    using Graph = StaticGraph<VertexAttributes, EdgeAttributes>;
    std::ifstream osmFile(osmFilename, std::ios::binary);
    std::ifstream refFile(refFilename, std::ios::binary);
    if (!osmFile.good())
      throw std::invalid_argument("file not found -- '" + osmFilename + "'");
    if (!refFile.good())
      throw std::invalid_argument("file not found -- '" + refFilename + "'");
    Graph osmGraph(osmFile);
    Graph refGraph(refFile);
    osmFile.close();
    refFile.close();

    // Build vp-trees for both graphs.
    std::vector<float> osmLat(osmGraph.numVertices());
    std::vector<float> osmLng(osmGraph.numVertices());
    std::vector<float> refLat(refGraph.numVertices());
    std::vector<float> refLng(refGraph.numVertices());
    FORALL_VERTICES(osmGraph, v) {
      osmLat[v] = osmGraph.latLng(v).latInDeg();
      osmLng[v] = osmGraph.latLng(v).lngInDeg();
    }
    FORALL_VERTICES(refGraph, v) {
      refLat[v] = refGraph.latLng(v).latInDeg();
      refLng[v] = refGraph.latLng(v).lngInDeg();
    }
    VPTree osmTree(osmLat, osmLng);
    VPTree refTree(refLat, refLng);

    // Match the vertices in both graphs with each other.
    int numMatchedVertices = 0;
    std::vector<int> osmToRefVertex(osmGraph.numVertices(), INVALID_VERTEX);
    FORALL_VERTICES(osmGraph, u) {
      const int v = findNearestNeighbor(refTree, osmGraph.latLng(u), maxDist);
      if (v != INVALID_INDEX) {
        const int w = findNearestNeighbor(osmTree, refGraph.latLng(v), maxDist);
        if (w == u) {
          ++numMatchedVertices;
          osmToRefVertex[u] = v;
        }
      }
    }

    // Open the output file.
    std::ofstream out(outfile + ".csv");
    if (!out.good())
      throw std::invalid_argument("file cannot be opened -- '" + outfile + ".csv'");
    out << "# OSM graph: " << osmFilename << std::endl;
    out << "# Ref graph: " << refFilename << std::endl;
    out << "# Maximum distance: " << maxDist << 'm' << std::endl;
    out << "osm_road_category,osm_speed_limit,free_flow_speed,capacity" << std::endl;

    // Match the edges in both graphs with each other and output all matched edges.
    int numMatchedEdges = 0;
    FORALL_VALID_EDGES(osmGraph, u1, e1) {
      const int v1 = osmGraph.edgeHead(e1);
      const int u2 = osmToRefVertex[u1];
      const int v2 = osmToRefVertex[v1];
      if (u2 == INVALID_VERTEX || v2 == INVALID_VERTEX)
        continue;
      if (osmGraph.uniqueEdgeBetween(u1, v1) == INVALID_EDGE)
        continue;
      const int e2 = refGraph.uniqueEdgeBetween(u2, v2);
      if (e2 == INVALID_EDGE)
        continue;

      out << osmGraph.osmRoadCategory(e1) << ',';
      if (osmGraph.speedLimit(e1) != -1)
        out << osmGraph.speedLimit(e1) << ',';
      else
        out << "NA,";
      out << refGraph.freeFlowSpeed(e2) << ',';
      out << refGraph.capacity(e2) / refGraph.numLanes(e2) << std::endl;
      ++numMatchedEdges;
    }

    const double matchedV1 = 100.0 * numMatchedVertices / osmGraph.numVertices();
    const double matchedV2 = 100.0 * numMatchedVertices / refGraph.numVertices();
    const double matchedE1 = 100.0 * numMatchedEdges / osmGraph.numEdges();
    const double matchedE2 = 100.0 * numMatchedEdges / refGraph.numEdges();
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Matched vertices in OSM graph: " << std::setw(6) << matchedV1 << '%' << '\n';
    std::cout << "Matched vertices in ref graph: " << std::setw(6) << matchedV2 << '%' << '\n';
    std::cout << "Matched edges in OSM graph: " << std::setw(6) << matchedE1 << '%' << '\n';
    std::cout << "Matched edges in ref graph: " << std::setw(6) << matchedE2 << '%' << '\n';
  } catch (std::exception& e) {
    std::cerr << argv[0] << ": " << e.what() << std::endl;
    std::cerr << "Try '" << argv[0] <<" -help' for more information." << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
