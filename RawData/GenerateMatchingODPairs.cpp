#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/dynamic_bitset.hpp>
#include <routingkit/geo_position_to_node.h>

#include "DataStructures/Geometry/Area.h"
#include "DataStructures/Geometry/LatLng.h"
#include "DataStructures/Geometry/Point.h"
#include "DataStructures/Graph/Attributes/LatLngAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "Tools/CommandLine/CommandLineParser.h"
#include "Tools/Constants.h"

inline void printUsage() {
  std::cout <<
      "Usage: GenerateMatchingODPairs -g1 <file> -g2 <file> -o1 <file> -o2 <file>\n"
      "This program takes two graphs as input, matches the vertices in both graphs with\n"
      "each other and outputs for each graph an OD-file representing the same journeys.\n"
      "  -n <num>          the number of OD-pairs to be generated\n"
      "  -d <meters>       only match vertices within a maximum distance of <meters>\n"
      "  -s <seed>         the seed for the random number generator\n"
      "  -g1 <file>        the first input graph in binary format\n"
      "  -g2 <file>        the second input graph in binary format\n"
      "  -p <file>         restrict the OD-pairs to an area from an OSM POLY file\n"
      "  -o1 <file>        place the output for the first graph in <file>\n"
      "  -o2 <file>        place the output for the second graph in <file>\n"
      "  -help             display this help and exit\n";
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

    const auto infile1 = clp.getValue<std::string>("g1");
    const auto infile2 = clp.getValue<std::string>("g2");
    const auto polyfile = clp.getValue<std::string>("p");
    const auto outfile1 = clp.getValue<std::string>("o1");
    const auto outfile2 = clp.getValue<std::string>("o2");
    const int numPairs = clp.getValue<int>("n");
    const int maxDist = clp.getValue<int>("d", 10);

    // Read both input graphs from disk.
    using Graph = StaticGraph<VertexAttrs<LatLngAttribute>>;
    std::ifstream in1(infile1, std::ios::binary);
    std::ifstream in2(infile2, std::ios::binary);
    if (!in1.good())
      throw std::invalid_argument("file not found -- '" + infile1 + "'");
    if (!in2.good())
      throw std::invalid_argument("file not found -- '" + infile2 + "'");
    Graph graph1(in1);
    Graph graph2(in2);
    in1.close();
    in2.close();

    // Build vp-trees for both graphs.
    std::vector<float> lat1(graph1.numVertices());
    std::vector<float> lat2(graph2.numVertices());
    std::vector<float> lng1(graph1.numVertices());
    std::vector<float> lng2(graph2.numVertices());
    FORALL_VERTICES(graph1, v) {
      lat1[v] = graph1.latLng(v).latInDeg();
      lng1[v] = graph1.latLng(v).lngInDeg();
    }
    FORALL_VERTICES(graph2, v) {
      lat2[v] = graph2.latLng(v).latInDeg();
      lng2[v] = graph2.latLng(v).lngInDeg();
    }
    VPTree tree1(lat1, lng1);
    VPTree tree2(lat2, lng2);

    // Check which vertices are inside the area of interest.
    boost::dynamic_bitset<> isVertexInsideArea1(graph1.numVertices());
    boost::dynamic_bitset<> isVertexInsideArea2(graph2.numVertices());
    isVertexInsideArea1.set();
    isVertexInsideArea2.set();
    if (!polyfile.empty()) {
      Area area;
      area.importFromOsmPolyFile(polyfile);
      const auto box = area.boundingBox();
      FORALL_VERTICES(graph1, v) {
        const Point p = {graph1.latLng(v).longitude(), graph1.latLng(v).latitude()};
        isVertexInsideArea1[v] = box.contains(p) && area.contains(p);
      }
      FORALL_VERTICES(graph2, v) {
        const Point p = {graph2.latLng(v).longitude(), graph2.latLng(v).latitude()};
        isVertexInsideArea2[v] = box.contains(p) && area.contains(p);
      }
    }

    // Match the vertices in both graphs with each other.
    std::vector<std::pair<int, int>> matchedVertices;
    FORALL_VERTICES(graph1, u) {
      if (isVertexInsideArea1[u]) {
        const int v = findNearestNeighbor(tree2, graph1.latLng(u), maxDist);
        if (v != INVALID_INDEX && isVertexInsideArea2[v]) {
          const int w = findNearestNeighbor(tree1, graph2.latLng(v), maxDist);
          if (w == u)
            matchedVertices.emplace_back(u, v);
        }
      }
    }

    std::cout << "Number of vertices in graph 1: " << graph1.numVertices() << std::endl;
    std::cout << "Number of vertices in graph 2: " << graph2.numVertices() << std::endl;
    std::cout << "Number of edges in graph 1: " << graph1.numEdges() << std::endl;
    std::cout << "Number of edges in graph 2: " << graph2.numEdges() << std::endl;
    const int numVerticesInsideArea1 = isVertexInsideArea1.count();
    const int numVerticesInsideArea2 = isVertexInsideArea2.count();
    const double matched1 = 100.0 * matchedVertices.size() / numVerticesInsideArea1;
    const double matched2 = 100.0 * matchedVertices.size() / numVerticesInsideArea2;
    if (!polyfile.empty()) {
      std::cout << "Vertices of interest in graph 1: " << numVerticesInsideArea1 << std::endl;
      std::cout << "Vertices of interest in graph 2: " << numVerticesInsideArea2 << std::endl;
    }
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Matched vertices in graph 1: " << std::setw(6) << matched1 << "%" << std::endl;
    std::cout << "Matched vertices in graph 2: " << std::setw(6) << matched2 << "%" << std::endl;

    // Open both output OD-files.
    std::ofstream out1(outfile1 + ".csv");
    std::ofstream out2(outfile2 + ".csv");
    if (!out1.good())
      throw std::invalid_argument("file cannot be opened -- '" + outfile1 + ".csv'");
    if (!out2.good())
      throw std::invalid_argument("file cannot be opened -- '" + outfile2 + ".csv'");
    out1 << "# Input graph: " << infile1 << std::endl;
    out2 << "# Input graph: " << infile2 << std::endl;
    out1 << "# Methodology: matched" << std::endl;
    out2 << "# Methodology: matched" << std::endl;
    out1 << "origin,destination" << std::endl;
    out2 << "origin,destination" << std::endl;

    // Pick the OD-pairs from the matched vertices uniformly at random.
    std::default_random_engine rand(clp.getValue<int>("s", 19900325));
    std::uniform_int_distribution<> dist(0, matchedVertices.size() - 1);
    for (int i = 0; i < numPairs; ++i) {
      const int o = dist(rand);
      const int d = dist(rand);
      out1 << matchedVertices[o].first << ',' << matchedVertices[d].first << std::endl;
      out2 << matchedVertices[o].second << ',' << matchedVertices[d].second << std::endl;
    }

    out1.close();
    out2.close();
  } catch (std::exception& e) {
    std::cerr << argv[0] << ": " << e.what() << std::endl;
    std::cerr << "Try '" << argv[0] <<" -help' for more information." << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
