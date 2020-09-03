#pragma once

#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <routingkit/osm_graph_builder.h>
#include <routingkit/tag_map.h>

#include "DataStructures/Containers/BitVector.h"
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

      const auto oneway = tags["oneway"];
      if (oneway) {
        if (stringEq(oneway, "no") || stringEq(oneway, "false") || stringEq(oneway, "0"))
          return RoadDirection::OPEN_IN_BOTH;
        else if (stringEq(oneway, "yes") || stringEq(oneway, "true") || stringEq(oneway, "1"))
          return RoadDirection::FORWARD;
        else if (stringEq(oneway, "-1") || stringEq(oneway, "reverse"))
          return RoadDirection::REVERSE;
      }

      const auto junction = tags["junction"];
      if (junction && stringEq(junction, "roundabout"))
        return RoadDirection::FORWARD;

      return roadDefaults.at(cat).direction;
    };

    // Returns the speed limit.
    auto getSpeedLimit = [&](const OsmRoadCategory cat, const RoutingKit::TagMap& tags) {
      assert(cat != OsmRoadCategory::ROAD);
      const auto maxspeed = tags["maxspeed"];
      if (maxspeed) {
        try {
          if (endsWith(maxspeed, "mph")) {
            std::string maxspeedWithoutUnit(maxspeed, std::strlen(maxspeed) - 3);
            trim(maxspeedWithoutUnit);
            return std::make_pair(lexicalCast<double>(maxspeedWithoutUnit) * 1.609344, true);
          } else {
            return std::make_pair(lexicalCast<double>(maxspeed), true);
          }
        } catch (std::logic_error& /*e*/) {}
      }
      return std::make_pair(static_cast<double>(roadDefaults.at(cat).speedLimit), false);
    };

    // Returns the number of lanes in the forward and reverse direction.
    auto getNumLanes =
        [&](const OsmRoadCategory cat, const RoadDirection dir, const RoutingKit::TagMap& tags) {
      assert(cat != OsmRoadCategory::ROAD);
      assert(dir != RoadDirection::CLOSED);
      const auto& defaults = roadDefaults.at(cat);
      auto numForwardLanes = defaults.numLanesPerDirection;
      auto numReverseLanes = defaults.numLanesPerDirection;

      const auto lanes = tags["lanes"];
      if (lanes) {
        try {
          auto numLanes = lexicalCast<double>(lanes);
          if (dir == RoadDirection::OPEN_IN_BOTH)
            numLanes /= 2;
          numForwardLanes = numLanes;
          numReverseLanes = numLanes;
        } catch (std::logic_error& /*e*/) {}
      }

      const auto lanesForward = tags["lanes:forward"];
      if (lanesForward) {
        try {
          auto numLanes = lexicalCast<double>(lanesForward);
          if (numLanes > 0)
            numForwardLanes = numLanes;
        } catch (std::logic_error& /*e*/) {}
      }

      const auto lanesBackward = tags["lanes:backward"];
      if (lanesBackward) {
        try {
          auto numLanes = lexicalCast<double>(lanesBackward);
          if (numLanes > 0)
            numReverseLanes = numLanes;
        } catch (std::logic_error& /*e*/) {}
      }

      return std::make_pair(numForwardLanes, numReverseLanes);
    };

    EnumParser<OsmRoadCategory> parseOsmRoadCategory;

    // Returns true if the specified way is open for vehicles.
    auto isWayOpenForVehicles = [&](uint64_t, const RoutingKit::TagMap& tags) {
      const auto highway = tags["highway"];
      if (highway)
        try {
          return parseOsmRoadCategory(highway) != OsmRoadCategory::ROAD;
        } catch (std::invalid_argument& /*e*/) {}
      return false;
    };

    // Invoked when a way is discovered that is open for vehicles.
    auto wayCallback = [&](uint64_t, const unsigned int seqId, const RoutingKit::TagMap& tags) {
      assert(seqId < wayCategory.size());
      assert(tags["highway"]);
      const auto cat = parseOsmRoadCategory(tags["highway"]);
      const auto dir = getRoadDirection(cat, tags);
      assert(dir != RoadDirection::CLOSED);

      wayCategory[seqId] = cat;
      auto bit = isWaySpeedLimitGiven[seqId];
      std::tie(waySpeedLimit[seqId], bit) = getSpeedLimit(cat, tags);
      std::tie(numForwardLanes[seqId], numReverseLanes[seqId]) = getNumLanes(cat, dir, tags);

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

    // Actually perform the work. Extract a graph from OSM data.
    const auto idMap = RoutingKit::load_osm_id_mapping_from_pbf(
        filename + ".osm.pbf", nullptr, isWayOpenForVehicles);
    const auto numWaysOpenForVehicles = idMap.is_routing_way.population_count();
    wayCategory.resize(numWaysOpenForVehicles);
    waySpeedLimit.resize(numWaysOpenForVehicles);
    isWaySpeedLimitGiven.resize(numWaysOpenForVehicles);
    numForwardLanes.resize(numWaysOpenForVehicles);
    numReverseLanes.resize(numWaysOpenForVehicles);
    osmGraph = load_osm_routing_graph_from_pbf(
        filename + ".osm.pbf", idMap, wayCallback, nullptr,
        nullptr, false, RoutingKit::OSMRoadGeometry::uncompressed);
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
    int speedLimit;           // The speed limit in km/h.
    int numLanesPerDirection; // The number of lanes per direction.
    int laneCapacity;         // The capacity per lane in veh/h.
    RoadDirection direction;  // The direction in which the road segment is open.
  };

  // A map that contains default road properties for a set of OSM road categories.
  inline static const std::unordered_map<OsmRoadCategory, RoadDefaults> roadDefaults = {
    {OsmRoadCategory::MOTORWAY,       {120, 2, 2000, RoadDirection::FORWARD     }},
    {OsmRoadCategory::MOTORWAY_LINK,  { 80, 1, 1500, RoadDirection::FORWARD     }},
    {OsmRoadCategory::TRUNK,          { 80, 1, 2000, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::TRUNK_LINK,     { 50, 1, 1500, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::PRIMARY,        { 80, 1, 1500, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::PRIMARY_LINK,   { 60, 1, 1500, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::SECONDARY,      { 30, 1,  800, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::SECONDARY_LINK, { 30, 1,  800, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::TERTIARY,       { 25, 1,  600, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::TERTIARY_LINK,  { 25, 1,  600, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::UNCLASSIFIED,   { 15, 1,  600, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::RESIDENTIAL,    { 15, 1,  600, RoadDirection::OPEN_IN_BOTH}},
    {OsmRoadCategory::LIVING_STREET,  { 10, 1,  300, RoadDirection::OPEN_IN_BOTH}},
  };

  RoutingKit::OSMRoutingGraph osmGraph;     // The graph extracted from OSM data.
  std::vector<OsmRoadCategory> wayCategory; // The OSM road category for each way.
  std::vector<double> waySpeedLimit;        // The speed limit for each way.
  BitVector isWaySpeedLimitGiven;           // Indicates for each way whether speed limit is given.
  std::vector<double> numForwardLanes;      // The number of forward lanes for each way.
  std::vector<double> numReverseLanes;      // The number of reverse lanes for each way.

  int currentVertex = -1; // The index of the current vertex in the OSM graph.
  int currentEdge = -1;   // The index of the current edge in the OSM graph.
  int currentTail = -1;   // The tail of the current edge.
};

// Returns the value of the LatLng attribute for the current vertex.
template <>
inline LatLngAttribute::Type OsmImporter::getValue<LatLngAttribute>() const {
  assert(currentVertex >= 0); assert(currentVertex < osmGraph.node_count());
  return {osmGraph.latitude[currentVertex], osmGraph.longitude[currentVertex]};
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
  const auto way = osmGraph.way[currentEdge];
  return osmGraph.is_arc_antiparallel_to_way[currentEdge]
      ? numReverseLanes[way] : numForwardLanes[way];
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
  const auto first = osmGraph.first_modelling_node[currentEdge];
  const auto last = osmGraph.first_modelling_node[currentEdge + 1];
  std::vector<LatLng> path(last - first);
  for (int v = first, i = 0; v < last; ++v, ++i)
    path[i] = {osmGraph.modelling_node_latitude[v], osmGraph.modelling_node_longitude[v]};
  return path;
}

// Returns the value of the speed limit attribute for the current edge.
template <>
inline SpeedLimitAttribute::Type OsmImporter::getValue<SpeedLimitAttribute>() const {
  assert(currentEdge >= 0); assert(currentEdge < osmGraph.arc_count());
  return std::round(waySpeedLimit[osmGraph.way[currentEdge]]);
}

// Returns the value of the capacity attribute for the current edge.
template <>
inline CapacityAttribute::Type OsmImporter::getValue<CapacityAttribute>() const {
  assert(currentEdge >= 0); assert(currentEdge < osmGraph.arc_count());
  auto laneCapacity = roadDefaults.at(wayCategory[osmGraph.way[currentEdge]]).laneCapacity;
  if (getValue<LengthAttribute>() < 100)
    laneCapacity *= 2;
  return std::round(getValue<NumLanesAttribute>() * laneCapacity);
}

// Returns the value of the free-flow speed attribute for the current edge.
template <>
inline FreeFlowSpeedAttribute::Type OsmImporter::getValue<FreeFlowSpeedAttribute>() const {
  assert(currentEdge >= 0); assert(currentEdge < osmGraph.arc_count());
  const auto way = osmGraph.way[currentEdge];
  if (isWaySpeedLimitGiven[way]) {
    return std::round(waySpeedLimit[way] < 51 ? 0.5 * waySpeedLimit[way] : waySpeedLimit[way]);
  } else {
    if (getValue<OsmRoadCategoryAttribute>() >= OsmRoadCategory::TRUNK &&
        getValue<OsmRoadCategoryAttribute>() <= OsmRoadCategory::TERTIARY_LINK &&
        getValue<LengthAttribute>() < 300) {
      return std::round((waySpeedLimit[way] - 10) / 300 * getValue<LengthAttribute>() + 10);
    }
    return std::round(waySpeedLimit[way]);
  }
}

// Returns the value of the travel time attribute for the current edge.
template <>
inline TravelTimeAttribute::Type OsmImporter::getValue<TravelTimeAttribute>() const {
  if (getValue<FreeFlowSpeedAttribute>() == 0)
    return INFTY;
  return std::round(36.0 * getValue<LengthAttribute>() / getValue<FreeFlowSpeedAttribute>());
}
