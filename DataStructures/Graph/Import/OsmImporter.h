#pragma once

#include <cassert>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <routingkit/osm_graph_builder.h>
#include <routingkit/tag_map.h>

#include "DataStructures/Graph/Attributes/CapacityAttribute.h"
#include "DataStructures/Graph/Attributes/FreeFlowSpeedAttribute.h"
#include "DataStructures/Graph/Attributes/LatLngAttribute.h"
#include "DataStructures/Graph/Attributes/LengthAttribute.h"
#include "DataStructures/Graph/Attributes/NumLanesAttribute.h"
#include "DataStructures/Graph/Attributes/OsmRoadCategoryAttribute.h"
#include "DataStructures/Graph/Attributes/RoadGeometryAttribute.h"
#include "DataStructures/Graph/Attributes/SpeedLimitAttribute.h"
#include "DataStructures/Graph/Attributes/TravelTimeAttribute.h"
#include "Tools/Constants.h"
#include "Tools/EnumParser.h"
#include "Tools/LexicalCast.h"
#include "Tools/StringHelpers.h"

// An importer for reading graphs from OpenStreetMap data. The OSM data must be given as a file in
// PBF format. Edge attributes depending on the mode of transportation (such as free-flow speed or
// road capacity) refer to cars. First, the class Graph repeatedly calls nextVertex to read the
// next vertex and fetches various vertex attributes. Then it repeatedly calls nextEdge to read
// the next edge and fetches various edge attributes.
class OsmImporter {
 public:
  // Opens the input file(s) and reads the header line(s).
  void init(const std::string& filename) {
    // Returns the direction in which the road segment is open.
    auto getRoadDirection = [&](const OsmRoadCategory cat, const RoutingKit::TagMap& tags) {
      assert(cat != OsmRoadCategory::ROAD);
      const char* oneway = tags["oneway"];
      if (oneway) {
        if (stringEq(oneway, "yes") || stringEq(oneway, "true") || stringEq(oneway, "1"))
          return RoadDirection::FORWARD;
        else if (stringEq(oneway, "-1") || stringEq(oneway, "reverse"))
          return RoadDirection::REVERSE;
        else if (stringEq(oneway, "no") || stringEq(oneway, "false") || stringEq(oneway, "0"))
          return RoadDirection::OPEN_IN_BOTH;
      }

      const char* junction = tags["junction"];
      if (junction && stringEq(junction, "roundabout"))
        return RoadDirection::FORWARD;
      return roadDefaults.at(cat).direction;
    };

    // Returns the speed limit.
    auto getSpeedLimit = [&](const OsmRoadCategory cat, const RoutingKit::TagMap& tags) {
      assert(cat != OsmRoadCategory::ROAD);
      const char* maxspeed = tags["maxspeed"];
      if (maxspeed) {
        try {
          if (endsWith(maxspeed, "mph")) {
            std::string maxspeedWithoutUnit(maxspeed, std::strlen(maxspeed) - 3);
            trim(maxspeedWithoutUnit);
            const double speedLimit = lexicalCast<double>(maxspeedWithoutUnit);
            if (speedLimit > 0)
              return speedLimit * 1.609344; // mph to km/h
          } else {
            const double speedLimit = lexicalCast<double>(maxspeed);
            if (speedLimit > 0)
              return speedLimit;
          }
        } catch (std::logic_error& /*e*/) {}
      }
      return static_cast<double>(roadDefaults.at(cat).speedLimit);
    };

    // Returns the number of lanes in the forward and reverse direction.
    auto getNumLanes =
        [&](const OsmRoadCategory cat, const RoadDirection dir, const RoutingKit::TagMap& tags) {
      assert(cat != OsmRoadCategory::ROAD);
      assert(dir != RoadDirection::CLOSED);
      const auto& defaults = roadDefaults.at(cat);
      double numLanesInForward = defaults.numLanesOnOneWay;
      double numLanesInReverse = defaults.numLanesOnOneWay;
      if (dir == RoadDirection::OPEN_IN_BOTH) {
        numLanesInForward = defaults.numLanesOnTwoWay;
        numLanesInReverse = defaults.numLanesOnTwoWay;
      }

      const char* lanes = tags["lanes"];
      if (lanes)
        try {
          const double totalNumLanes = lexicalCast<double>(lanes);
          if (totalNumLanes > 0) {
            // Get the number of lanes in the forward direction.
            const char* lanesForward = tags["lanes:forward"];
            if (lanesForward) {
              const double parsedNumLanesInForward = lexicalCast<double>(lanesForward);
              if (parsedNumLanesInForward > 0)
                numLanesInForward = parsedNumLanesInForward;
            } else {
              // If this is a two-way segment, distribute the lanes evenly over both directions.
              if (dir != RoadDirection::OPEN_IN_BOTH)
                numLanesInForward = totalNumLanes;
              else
                numLanesInForward = totalNumLanes / 2;
            }

            // Get the number of lanes in the reverse direction.
            const char* lanesBackward = tags["lanes:backward"];
            if (lanesBackward) {
              const double parsedNumLanesInReverse = lexicalCast<double>(lanesBackward);
              if (parsedNumLanesInReverse > 0)
                numLanesInReverse = parsedNumLanesInReverse;
            } else {
              // If this is a two-way segment, distribute the lanes evenly over both directions.
              if (dir != RoadDirection::OPEN_IN_BOTH)
                numLanesInReverse = totalNumLanes;
              else
                numLanesInReverse = totalNumLanes / 2;
            }
          }
        } catch (std::logic_error& /*e*/) {}
      return std::make_pair(numLanesInForward, numLanesInReverse);
    };

    EnumParser<OsmRoadCategory> parseOsmRoadCategory;

    // Returns true if the OSM way is open for cars.
    auto isWayOpenForCars = [&](const uint64_t /*origId*/, const RoutingKit::TagMap& tags) {
      const char* access = tags["access"];
      if (access && stringEq(access, "no"))
        return false;

      const char* highway = tags["highway"];
      if (highway)
        try {
          return parseOsmRoadCategory(highway) != OsmRoadCategory::ROAD;
        } catch (std::invalid_argument& /*e*/) {}
      return false;
    };

    // Invoked when a way is discovered that is open for cars.
    auto wayCallback = [=, &parseOsmRoadCategory]
        (const uint64_t /*origId*/, const unsigned int seqId, const RoutingKit::TagMap& tags)
        -> RoutingKit::OSMWayDirectionCategory {
      assert(seqId < wayCategory.size());
      assert(tags["highway"]);
      const auto cat = parseOsmRoadCategory(tags["highway"]);
      const auto dir = getRoadDirection(cat, tags);
      assert(dir != RoadDirection::CLOSED);

      wayCategory[seqId] = cat;
      waySpeed[seqId] = getSpeedLimit(cat, tags);
      std::tie(numLanesInForward[seqId], numLanesInReverse[seqId]) = getNumLanes(cat, dir, tags);

      switch (dir) {
        case RoadDirection::OPEN_IN_BOTH:
          return RoutingKit::OSMWayDirectionCategory::open_in_both;
        case RoadDirection::FORWARD:
          return RoutingKit::OSMWayDirectionCategory::only_open_forwards;
        case RoadDirection::REVERSE:
          return RoutingKit::OSMWayDirectionCategory::only_open_backwards;
        default:
          assert(false);
          return RoutingKit::OSMWayDirectionCategory::closed;
      }
    };

    // Eventually, actually perform the work. Extract a graph from OSM data.
    const auto mapping =
        RoutingKit::load_osm_id_mapping_from_pbf(filename + ".osm.pbf", nullptr, isWayOpenForCars);
    const int numWaysOpenForCars = mapping.is_routing_way.population_count();
    wayCategory.resize(numWaysOpenForCars);
    waySpeed.resize(numWaysOpenForCars);
    numLanesInForward.resize(numWaysOpenForCars);
    numLanesInReverse.resize(numWaysOpenForCars);
    osmGraph = load_osm_routing_graph_from_pbf(
        filename + ".osm.pbf", mapping, wayCallback, nullptr, nullptr,
        false, RoutingKit::OSMRoadGeometry::uncompressed);
  }

  // Returns the number of vertices in the graph, or 0 if the number is not yet known.
  int numVertices() const {
    return osmGraph.node_count();
  }

  // Returns the number of edges in the graph, or 0 if the number is not yet known.
  int numEdges() const {
    return osmGraph.arc_count();
  }

  // Reads the next vertex from disk. Returns false if there are no more vertices.
  bool nextVertex() {
    return ++currentVertex < osmGraph.node_count();
  }

  // Returns the ID of the current vertex. Vertices must have sequential IDs from 0 to n − 1.
  int vertexId() const {
    return currentVertex;
  }

  // Reads the next edge from disk. Returns false if there are no more edges.
  bool nextEdge() {
    ++currentEdge;
    while (currentTail < static_cast<int>(osmGraph.node_count()) &&
        osmGraph.first_out[currentTail + 1] == currentEdge)
      ++currentTail;
    return currentEdge < osmGraph.arc_count();
  }

  // Returns the tail vertex of the current edge.
  int edgeTail() const {
    return currentTail;
  }

  // Returns the head vertex of the current edge.
  int edgeHead() const {
    assert(currentEdge >= 0); assert(currentEdge < osmGraph.arc_count());
    return osmGraph.head[currentEdge];
  }

  // Returns the value of the specified attribute for the current vertex/edge, or the attribute's
  // default value if the attribute is not part of the file format.
  template <typename Attr>
  typename Attr::Type getValue() const {
    return Attr::defaultValue();
  }

  // Closes the input file(s).
  void close() {}

 private:
  // A struct that carries standard values for a set of static road properties.
  struct RoadDefaults {
    int speedLimit;          // The speed limit in km/h.
    double freeFlowFactor;   // free-flow speed = free-flow factor * speed limit
    int numLanesOnOneWay;    // The number of lanes on one-way road segments.
    int numLanesOnTwoWay;    // The number of lanes in each direction on two-way road segments.
    int laneCapacity;        // The capacity per lane in veh/h.
    RoadDirection direction; // The direction in which the road segment is open.
  };

  // A map that carries default road properties for a set of OSM road categories.
  const std::unordered_map<OsmRoadCategory, RoadDefaults> roadDefaults = {
#ifndef CALIBRATE_OSM
    {OsmRoadCategory::MOTORWAY,       {120, 1.0, 2, 2, 800, RoadDirection::FORWARD}},
    {OsmRoadCategory::TRUNK,          {80,  1.0, 2, 1, 800, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::PRIMARY,        {80,  1.0, 2, 1, 600, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::SECONDARY,      {60,  1.0, 2, 1, 400, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::TERTIARY,       {45,  1.0, 1, 1, 240, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::UNCLASSIFIED,   {45,  1.0, 1, 1, 240, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::RESIDENTIAL,    {30,  1.0, 1, 1, 240, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::MOTORWAY_LINK,  {80,  1.0, 1, 1, 600, RoadDirection::FORWARD}},
    {OsmRoadCategory::TRUNK_LINK,     {50,  1.0, 1, 1, 600, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::PRIMARY_LINK,   {60,  1.0, 1, 1, 600, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::SECONDARY_LINK, {60,  1.0, 1, 1, 400, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::TERTIARY_LINK,  {45,  1.0, 1, 1, 240, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::LIVING_STREET,  {15,  1.0, 1, 1, 120, RoadDirection::OPEN_IN_BOTH}}
#else
    {OsmRoadCategory::MOTORWAY,       {-1, 1.0, 2, 2, -1, RoadDirection::FORWARD}},
    {OsmRoadCategory::TRUNK,          {-1, 1.0, 2, 1, -1, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::PRIMARY,        {-1, 1.0, 2, 1, -1, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::SECONDARY,      {-1, 1.0, 2, 1, -1, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::TERTIARY,       {-1, 1.0, 1, 1, -1, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::UNCLASSIFIED,   {-1, 1.0, 1, 1, -1, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::RESIDENTIAL,    {-1, 1.0, 1, 1, -1, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::MOTORWAY_LINK,  {-1, 1.0, 1, 1, -1, RoadDirection::FORWARD}},
    {OsmRoadCategory::TRUNK_LINK,     {-1, 1.0, 1, 1, -1, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::PRIMARY_LINK,   {-1, 1.0, 1, 1, -1, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::SECONDARY_LINK, {-1, 1.0, 1, 1, -1, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::TERTIARY_LINK,  {-1, 1.0, 1, 1, -1, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::LIVING_STREET,  {-1, 1.0, 1, 1, -1, RoadDirection::OPEN_IN_BOTH}}
#endif
  };

  RoutingKit::OSMRoutingGraph osmGraph;     // The graph extracted from OSM data.
  std::vector<OsmRoadCategory> wayCategory; // The OSM road category for each way.
  std::vector<double> waySpeed;             // The speed limit for each way.
  std::vector<double> numLanesInForward;    // The number of forward lanes for each way.
  std::vector<double> numLanesInReverse;    // The number of reverse lanes for each way.

  int currentVertex = -1; // The index of the current vertex in the OSM graph.
  int currentTail = -1;   // The tail vertex of the current edge.
  int currentEdge = -1;   // The index of the current edge in the OSM graph.
};

// Returns the value of the LatLng attribute for the current vertex.
template <>
inline LatLngAttribute::Type OsmImporter::getValue<LatLngAttribute>() const {
  assert(currentVertex >= 0); assert(currentVertex < osmGraph.node_count());
  return {osmGraph.latitude[currentVertex], osmGraph.longitude[currentVertex]};
}

// Returns the value of the capacity attribute for the current edge.
template <>
inline CapacityAttribute::Type OsmImporter::getValue<CapacityAttribute>() const {
  assert(currentEdge >= 0); assert(currentEdge < osmGraph.arc_count());
  const int way = osmGraph.way[currentEdge];
  const double numLanes = osmGraph.is_arc_antiparallel_to_way[currentEdge] ?
      numLanesInReverse[way] : numLanesInForward[way];
  const int laneCapacity = roadDefaults.at(wayCategory[way]).laneCapacity;
  return std::round(numLanes * laneCapacity);
}

// Returns the value of the free-flow speed attribute for the current edge.
template <>
inline FreeFlowSpeedAttribute::Type OsmImporter::getValue<FreeFlowSpeedAttribute>() const {
  assert(currentEdge >= 0); assert(currentEdge < osmGraph.arc_count());
  const int way = osmGraph.way[currentEdge];
  return std::round(roadDefaults.at(wayCategory[way]).freeFlowFactor * waySpeed[way]);
}

// Returns the value of the length attribute for the current edge.
template <>
inline LengthAttribute::Type OsmImporter::getValue<LengthAttribute>() const {
  assert(currentEdge >= 0); assert(currentEdge < osmGraph.arc_count());
  return osmGraph.geo_distance[currentEdge];
}

// Returns the value of the number of lanes attribute for the current edge.
template <>
inline NumLanesAttribute::Type OsmImporter::getValue<NumLanesAttribute>() const {
  assert(currentEdge >= 0); assert(currentEdge < osmGraph.arc_count());
  const int way = osmGraph.way[currentEdge];
  return osmGraph.is_arc_antiparallel_to_way[currentEdge] ?
      numLanesInReverse[way] : numLanesInForward[way];
}

// Returns the value of the OSM road category attribute for the current edge.
template <>
inline OsmRoadCategoryAttribute::Type OsmImporter::getValue<OsmRoadCategoryAttribute>() const {
  assert(currentEdge >= 0); assert(currentEdge < osmGraph.arc_count());
  return wayCategory[osmGraph.way[currentEdge]];
}

// Returns the value of the road geometry attribute for the current edge.
template <>
inline RoadGeometryAttribute::Type OsmImporter::getValue<RoadGeometryAttribute>() const {
  assert(currentEdge >= 0); assert(currentEdge < osmGraph.arc_count());
  const int first = osmGraph.first_modelling_node[currentEdge];
  const int last = osmGraph.first_modelling_node[currentEdge + 1];
  std::vector<LatLng> path(last - first);
  for (int v = first, i = 0; v != last; ++v, ++i)
    path[i] = {osmGraph.modelling_node_latitude[v], osmGraph.modelling_node_longitude[v]};
  return path;
}

// Returns the value of the speed limit attribute for the current edge.
template <>
inline SpeedLimitAttribute::Type OsmImporter::getValue<SpeedLimitAttribute>() const {
  assert(currentEdge >= 0); assert(currentEdge < osmGraph.arc_count());
  return std::round(waySpeed[osmGraph.way[currentEdge]]);
}

// Returns the value of the travel time attribute for the current edge.
template <>
inline TravelTimeAttribute::Type OsmImporter::getValue<TravelTimeAttribute>() const {
  assert(currentEdge >= 0); assert(currentEdge < osmGraph.arc_count());
  return std::round(36 * osmGraph.geo_distance[currentEdge] / waySpeed[osmGraph.way[currentEdge]]);
}
