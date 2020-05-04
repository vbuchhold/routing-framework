#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <csv.h>

#include "Algorithms/DemandCalculation/PopulationAssignment.h"
#include "Algorithms/Dijkstra/Dijkstra.h"
#include "DataStructures/Containers/BitVector.h"
#include "DataStructures/Geometry/Area.h"
#include "DataStructures/Geometry/KDTree.h"
#include "DataStructures/Geometry/Point.h"
#include "DataStructures/Graph/Attributes/LatLngAttribute.h"
#include "DataStructures/Graph/Attributes/LengthAttribute.h"
#include "DataStructures/Graph/Attributes/NumOpportunitiesAttribute.h"
#include "DataStructures/Graph/Attributes/PopulationAttribute.h"
#include "DataStructures/Graph/Attributes/TravelTimeAttribute.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"
#include "Tools/CommandLine/CommandLineParser.h"
#include "Tools/CommandLine/ProgressBar.h"
#include "Tools/StringHelpers.h"
#include "Tools/Timer.h"

inline void printUsage() {
  std::cout <<
      "Usage: MeanDistancePerRank -f <fmt> -g <file> -pop <file> -o <file>\n"
      "Computes the expected shortest-path distance from a uniform random inhabitant to\n"
      "the opportunity with rank r, for each rank r. Opportunities are ranked according\n"
      "to either the shortest-path metric or the geometric metric.\n"
      "  -len              use physical lengths as metric (default: travel times)\n"
      "  -r <num>          use Moore neighborhoods of max range <num> (defaults to 1)\n"
      "  -d <meters>       do not assign POIs to vertices farther than <meters>\n"
      "  -f <fmt>          format of the population grid file\n"
      "                      possible values: DE EU\n"
      "  -g <file>         input graph in binary format\n"
      "  -a <file>         restrict origins and destinations to polygonal study area\n"
      "  -pop <file>       population grid\n"
      "  -poi <file>       use density of POIs as proxy for trip attraction rates\n"
      "  -o <file>         place output in <file>\n"
      "  -help             display this help and exit\n";
}

// A graph data structure encompassing all required attributes.
using VertexAttributes = VertexAttrs<
    LatLngAttribute, NumOpportunitiesAttribute, PopulationAttribute>;
using EdgeAttributes = EdgeAttrs<LengthAttribute, TravelTimeAttribute>;
using GraphT = StaticGraph<VertexAttributes, EdgeAttributes>;

// Writes the specified distances to file. The rankingMetric parameter should be "Dij" or "geo".
inline void writeDistPerRankToFile(
    const std::vector<int64_t>& distPerRank, const std::string& rankingMetric,
    const CommandLineParser& clp) {
  const auto maxRange = clp.getValue<int>("r", 1);
  const auto maxDistance = clp.getValue<int>("d", 200);
  const auto graphFileName = clp.getValue<std::string>("g");
  const auto areaFileName = clp.getValue<std::string>("a");
  const auto popFileName = clp.getValue<std::string>("pop");
  const auto poiFileName = clp.getValue<std::string>("poi");
  auto outputFileName = clp.getValue<std::string>("o");
  if (endsWith(outputFileName, ".csv"))
    outputFileName.erase(outputFileName.size() - 4);
  outputFileName += '.' + rankingMetric + ".csv";

  std::ofstream outputFile(outputFileName);
  if (!outputFile.good())
    throw std::invalid_argument("file cannot be opened -- '" + outputFileName + "'");

  outputFile << "# Graph: " << graphFileName << '\n';
  if (!areaFileName.empty())
    outputFile << "# Area: " << areaFileName << '\n';
  outputFile << "# POP: " << popFileName << " (r=" << maxRange << ")\n";
  if (!areaFileName.empty())
    outputFile << "# POI: " << poiFileName << " (d=" << maxDistance << ")\n";
  outputFile << rankingMetric << "_rank,expected_distance\n";

  for (auto i = 0; i < distPerRank.size(); ++i)
    outputFile << i << ',' << distPerRank[i] << '\n';
}

// Compute and output the expected shortest-path distances for Dijkstra and geometric ranks.
inline void computeAndOutputDistPerRank(const GraphT& graph, const CommandLineParser& clp) {
  Timer timer;
  std::cout << "Computing distances: ";
  ProgressBar bar(graph.numVertices());

  auto population = 0.0;
  auto numOpportunities = 0.0;
  FORALL_VERTICES(graph, v) {
    population += graph.population(v);
    numOpportunities += graph.numOpportunities(v);
  }

  std::vector<int64_t> distPerDijRank(numOpportunities);
  std::vector<int64_t> distPerGeoRank(numOpportunities);

  #pragma omp parallel
  {
    struct Vertex {
      int id;
      int dijDistance;
      int64_t geoDistance;
    };
    std::vector<Vertex> verticesOrderedByRank(graph.numVertices());

    std::vector<int64_t> localDistPerDijRank(numOpportunities);
    std::vector<int64_t> localDistPerGeoRank(numOpportunities);

    Dijkstra<GraphT, TravelTimeAttribute, BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>> dij(graph);

    #pragma omp for schedule(static, 1) nowait
    FORALL_VERTICES(graph, s) {
      ++bar;
      if (graph.population(s) == 0)
        continue;
      dij.run(s);

      const auto srcCoordinate = Point(graph.latLng(s).longitude(), graph.latLng(s).latitude());
      FORALL_VERTICES(graph, v) {
        const auto coordinate = Point(graph.latLng(v).longitude(), graph.latLng(v).latitude());
        auto& vertex = verticesOrderedByRank[v];
        vertex.id = v;
        vertex.dijDistance = dij.getDistance(v);
        vertex.geoDistance = srcCoordinate.getSquaredEuclideanDistanceTo(coordinate);
      }

      // Rank the vertices according to the shortest-path metric.
      std::sort(verticesOrderedByRank.begin(), verticesOrderedByRank.end(), [&](auto u, auto v) {
        return u.dijDistance < v.dijDistance;
      });

      for (auto i = 0, rank = 0; i < graph.numVertices(); ++i) {
        const auto v = verticesOrderedByRank[i];
        for (auto j = 0; j < graph.numOpportunities(v.id); ++j, ++rank) {
          assert(rank < numOpportunities);
          localDistPerDijRank[rank] += graph.population(s) * v.dijDistance;
        }
      }

      // Rank the vertices according to the geometric metric.
      std::sort(verticesOrderedByRank.begin(), verticesOrderedByRank.end(), [&](auto u, auto v) {
        return u.geoDistance < v.geoDistance || (u.geoDistance == v.geoDistance && u.id < v.id);
      });

      for (auto i = 0, rank = 0; i < graph.numVertices(); ++i) {
        const auto v = verticesOrderedByRank[i];
        for (auto j = 0; j < graph.numOpportunities(v.id); ++j, ++rank) {
          assert(rank < numOpportunities);
          localDistPerGeoRank[rank] += graph.population(s) * v.dijDistance;
        }
      }
    }

    #pragma omp critical (mergeResults)
    {
      for (auto i = 0; i < numOpportunities; ++i) {
        distPerDijRank[i] += localDistPerDijRank[i];
        distPerGeoRank[i] += localDistPerGeoRank[i];
      }
    }
  }

  for (auto i = 0; i < numOpportunities; ++i) {
    distPerDijRank[i] = std::round(distPerDijRank[i] / population);
    distPerGeoRank[i] = std::round(distPerGeoRank[i] / population);
  }

  bar.finish();
  std::cout << "done (" << timer.elapsed<std::chrono::seconds>() << "s).\n";

  std::cout << "Writing distances to file..." << std::flush;
  writeDistPerRankToFile(distPerDijRank, "Dij", clp);
  writeDistPerRankToFile(distPerGeoRank, "geo", clp);
  std::cout << " done.\n";
}

int main(int argc, char* argv[]) {
  try {
    CommandLineParser clp(argc, argv);
    if (clp.isSet("help")) {
      printUsage();
      return EXIT_SUCCESS;
    }

    // Parse the command-line options.
    const auto useLengths = clp.isSet("len");
    const auto maxRange = clp.getValue<int>("r", 1);
    const auto maxDistance = clp.getValue<int>("d", 200);
    const auto popFileFormat = clp.getValue<std::string>("f");
    const auto graphFileName = clp.getValue<std::string>("g");
    const auto areaFileName = clp.getValue<std::string>("a");
    const auto popFileName = clp.getValue<std::string>("pop");
    const auto poiFileName = clp.getValue<std::string>("poi");

    // Validate command-line options.
    if (maxRange < 0)
      throw std::invalid_argument("max range is smaller than 0");
    if (maxDistance < 0)
      throw std::invalid_argument("max distance is smaller than 0");

    // Read the graph from file.
    std::cout << "Reading graph from file..." << std::flush;
    std::ifstream graphFile(graphFileName, std::ios::binary);
    if (!graphFile.good())
      throw std::invalid_argument("file not found -- '" + graphFileName + "'");
    GraphT graph(graphFile);
    graphFile.close();
    if (useLengths)
      FORALL_EDGES(graph, e)
        graph.travelTime(e) = graph.length(e);
    std::cout << " done.\n";

    // Read the study area from OSM POLY file.
    BitVector isVertexInStudyArea(graph.numVertices(), true);
    if (!areaFileName.empty()) {
      std::cout << "Reading study area from OSM POLY file..." << std::flush;
      Area studyArea;
      studyArea.importFromOsmPolyFile(areaFileName);
      const auto box = studyArea.boundingBox();
      FORALL_VERTICES(graph, u) {
        const Point p(graph.latLng(u).longitude(), graph.latLng(u).latitude());
        isVertexInStudyArea[u] = box.contains(p) && studyArea.contains(p);
      }
      std::cout << " done.\n";
    }

    // Assign the population grid to the graph.
    if (popFileFormat == "DE") {
      using PopFileReaderT = io::CSVReader<2, io::trim_chars<>, io::no_quote_escape<';'>>;
      using PopAssignmentT = PopulationAssignment<GraphT, PopFileReaderT>;
      PopFileReaderT popFileReader(popFileName);
      popFileReader.read_header(io::ignore_extra_column, "Gitter_ID_100m", "Einwohner");
      PopAssignmentT assign(graph, isVertexInStudyArea, popFileReader, 100, maxRange);
      assign.run(true);
    } else if (popFileFormat == "EU") {
      using PopFileReaderT = io::CSVReader<2>;
      using PopAssignmentT = PopulationAssignment<GraphT, PopFileReaderT>;
      PopFileReaderT popFileReader(popFileName);
      popFileReader.read_header(io::ignore_extra_column, "GRD_ID", "TOT_P");
      PopAssignmentT assign(graph, isVertexInStudyArea, popFileReader, 1000, maxRange);
      assign.run(true);
    } else {
      throw std::invalid_argument("invalid population file format -- '" + popFileFormat + "'");
    }

    // Assign POIs to vertices (or use the density of population as a proxy for trip attraction).
    if (!poiFileName.empty()) {
      Timer timer;
      std::cout << "Assigning POIs to vertices..." << std::flush;
      int assignedPoi = 0;
      int unassignedPoi = 0;

      std::vector<int> verticesInStudyArea(isVertexInStudyArea.cardinality());
      for (auto i = 0, v = isVertexInStudyArea.firstSetBit(); i < verticesInStudyArea.size(); ++i) {
        verticesInStudyArea[i] = v;
        v = isVertexInStudyArea.nextSetBit(v);
      }

      std::vector<Point> points(verticesInStudyArea.size());
      for (auto i = 0; i < verticesInStudyArea.size(); ++i) {
        const auto& latLng = graph.latLng(verticesInStudyArea[i]);
        points[i] = {latLng.longitude(), latLng.latitude()};
      }
      KDTree tree(points);

      double lng, lat;
      using TrimPolicy = io::trim_chars<>;
      using QuotePolicy = io::no_quote_escape<','>;
      using OverflowPolicy = io::throw_on_overflow;
      using CommentPolicy = io::single_line_comment<'#'>;
      using PoiReader = io::CSVReader<2, TrimPolicy, QuotePolicy, OverflowPolicy, CommentPolicy>;
      PoiReader poiReader(poiFileName);
      poiReader.read_header(io::ignore_extra_column, "longitude", "latitude");
      while (poiReader.read_row(lng, lat)) {
        const LatLng query(lat, lng);
        const auto p = tree.findClosestPoint(Point(query.longitude(), query.latitude()));
        if (graph.latLng(verticesInStudyArea[p]).getGreatCircleDistanceTo(query) <= maxDistance) {
          ++graph.numOpportunities(verticesInStudyArea[p]);
          ++assignedPoi;
        } else {
          ++unassignedPoi;
        }
      }

      std::cout << " done (" << timer.elapsed() << "ms).\n";
      std::cout << "  POIs that were assigned: " << assignedPoi << "\n";
      std::cout << "  POIs that could not be assigned: " << unassignedPoi << "\n";
    } else {
      FORALL_VERTICES(graph, v)
        graph.numOpportunities(v) = graph.population(v);
    }

    computeAndOutputDistPerRank(graph, clp);
  } catch (std::exception& e) {
    std::cerr << argv[0] << ": " << e.what() << '\n';
    std::cerr << "Try '" << argv[0] <<" -help' for more information.\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
