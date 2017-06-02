#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>

#include "DataStructures/Geometry/LatLng.h"
#include "DataStructures/Graph/Attributes/CapacityAttribute.h"
#include "DataStructures/Graph/Attributes/LatLngAttribute.h"
#include "DataStructures/Graph/Attributes/LengthAttribute.h"
#include "DataStructures/Graph/Attributes/TravelTimeAttribute.h"
#include "DataStructures/Graph/Attributes/XatfRoadCategoryAttribute.h"
#include "Tools/LexicalCast.h"

// An importer for reading graph files in the XATF file format. First, the Graph class repeatedly
// calls nextVertex() to read the next vertex from disk and fetches various attributes of the
// vertex. Then, it repeatedly calls nextEdge() to read the next edge from disk and fetches various
// attributes of the edge.
class XatfImporter {
 public:
  // Constructs an XatfImporter.
  XatfImporter() : returnReverse(false), nextVertexId(0) {}

  // Constructs an XatfImporter with the same configuration as the specified XatfImporter.
  XatfImporter(const XatfImporter& /*other*/) : XatfImporter() {}

  // Opens the input file(s) and reads the header line(s).
  void init(const std::string& filename) {
    verticesFile.open(filename + ".nfx");
    assert(verticesFile);
    edgesFile.open(filename + ".sfx");
    assert(edgesFile);
  }

  // Returns the number of vertices in the graph, or 0 if the number is not yet known.
  int numVertices() const {
    return 0;
  }

  // Returns the number of edges in the graph, or 0 if the number is not yet known.
  int numEdges() const {
    return 0;
  }

  // Reads the next vertex from disk. Returns false if there are no more vertices, true otherwise.
  bool nextVertex() {
    // Read the next vertex from disk.
    int fieldsRead = nextRecord(verticesFile);
    if (fieldsRead == 0)
      return false;
    assert(fieldsRead == NUM_VERTEX_FIELDS);

    // Map the original vertex ID to a sequential ID.
    assert(origToNewIds.count(currentRecord[VERTEX_ID]) == 0);
    origToNewIds[currentRecord[VERTEX_ID]] = nextVertexId++;

    // Store the LatLng of the vertex (possibly needed to compute edge lengths).
    const int lat = lexicalCast<int>(currentRecord[LATITUDE]) * LATLNG_FACTOR;
    const int lng = lexicalCast<int>(currentRecord[LONGITUDE]) * LATLNG_FACTOR;
    latLngs.emplace_back(lat, lng);
    assert(latLngs.size() == nextVertexId);
    return true;
  }

  // Returns the ID of the current vertex. Vertices have to have sequential IDs from 0 to n − 1.
  int vertexId() const {
    return origToNewIds.at(currentRecord[VERTEX_ID]);
  }

  // Reads the next edge from disk. Returns false if there are no more edges, true otherwise.
  bool nextEdge() {
    if (returnReverse) {
      std::swap(currentRecord[EDGE_TAIL], currentRecord[EDGE_HEAD]);
      returnReverse = false;
      return true;
    }

    // Read the next edge from disk.
    int fieldsRead = nextRecord(edgesFile);
    if (fieldsRead == 0)
      return false;
    assert(fieldsRead == NUM_EDGE_FIELDS);

    currentEdge.length = lexicalCast<int>(currentRecord[LENGTH]);
    currentEdge.xatfRoadCategory = lexicalCast<XatfRoadCategory>(currentRecord[XATF_ROAD_CATEGORY]);

    // Validate edge record.
    assert(origToNewIds.count(currentRecord[EDGE_TAIL]) == 1);
    assert(origToNewIds.count(currentRecord[EDGE_HEAD]) == 1);
    assert(currentEdge.length >= 0);

    // Handle edge direction. We regard closed edges as bidirected to obtain a larger graph.
    Direction dir = static_cast<Direction>(lexicalCast<int>(currentRecord[DIRECTION]));
    assert(dir == BIDIRECTED || dir == FORWARD || dir == REVERSE || dir == CLOSED);
    returnReverse = dir == BIDIRECTED || dir == CLOSED;
    if (dir == REVERSE)
      std::swap(currentRecord[EDGE_TAIL], currentRecord[EDGE_HEAD]);
    return true;
  }

  // Returns the tail vertex of the current edge.
  int edgeTail() const {
    return origToNewIds.at(currentRecord[EDGE_TAIL]);
  }

  // Returns the head vertex of the current edge.
  int edgeHead() const {
    return origToNewIds.at(currentRecord[EDGE_HEAD]);
  }

  // Returns the value of the specified attribute for the current vertex/edge, or the attribute's
  // default value if the attribute is not part of the file format.
  template <typename Attr>
  typename Attr::Type getValue() const {
    return Attr::DEFAULT_VALUE;
  }

  // Closes the input file(s).
  void close() {
    verticesFile.close();
    edgesFile.close();
  }

 private:
  // An edge record in the XATF file format.
  struct EdgeRecord {
    int length;
    XatfRoadCategory xatfRoadCategory;
  };

  // These enums associate field names with indices in the vertex/edge record.
  enum VertexField {
    VERTEX_ID = 0,
    LONGITUDE = 3,
    LATITUDE  = 4,
  };

  enum EdgeField {
    EDGE_TAIL = 1,
    EDGE_HEAD = 2,
    LENGTH    = 4,
    DIRECTION = 6,
    XATF_ROAD_CATEGORY = 9,
  };

  // This enum specifies possible values for the "direction" edge field.
  enum Direction {
    BIDIRECTED = 0,
    FORWARD = 1,
    REVERSE = 2,
    CLOSED  = 3,
  };

  // The number of fields in each vertex/edge record.
  static constexpr int NUM_VERTEX_FIELDS = 5;
  static constexpr int NUM_EDGE_FIELDS   = 10;
  static constexpr int MAX_NUM_FIELDS    = std::max(NUM_VERTEX_FIELDS, NUM_EDGE_FIELDS);

  // Latitudes and longitudes are stored as integers in units of 1/LATLNG_PRECISION degrees.
  static constexpr int LATLNG_PRECISION = 100000;
  static constexpr int LATLNG_FACTOR = LatLng::PRECISION / LATLNG_PRECISION;

  // Travel times (of ferries) are stored as integers in units of TRAVEL_TIME_PRECISION seconds.
  static constexpr int TRAVEL_TIME_PRECISION = 60;
  static constexpr int TRAVEL_TIME_FACTOR = TRAVEL_TIME_PRECISION * 10;

  // Speed limits in km/h for the XATF road categories.
  static constexpr const int (&SPEED_LIMITS)[] = {
    -1,
    130, // Motorway, fast.
    120, // Motorway, medium.
    110, // Motorway, slow.
    100, // National road, fast.
    90,  // National road, medium.
    80,  // National road, slow.
    70,  // Regional road, fast.
    60,  // Regional road, medium.
    50,  // Regional road, slow.
    40,  // Urban street, fast.
    30,  // Urban street, medium.
    20,  // Urban street, slow.
    -1,  // Ferry.
    -1,  // Unused.
    10,  // Forest road.
  };

  // Numbers of lanes for the XATF road categories.
  static constexpr const int (&NUM_LANES)[] = {
    -1,
    2,  // Motorway, fast.
    2,  // Motorway, medium.
    2,  // Motorway, slow.
    2,  // National road, fast.
    1,  // National road, medium.
    1,  // National road, slow.
    1,  // Regional road, fast.
    1,  // Regional road, medium.
    1,  // Regional road, slow.
    1,  // Urban street, fast.
    1,  // Urban street, medium.
    1,  // Urban street, slow.
    -1, // Ferry.
    -1, // Unused.
    1,  // Forest road.
  };

  static constexpr int VEHICLE_LENGTH = 5; // The average vehicle length in meters.

  using RecordT = std::array<std::string, MAX_NUM_FIELDS>;
  using VertexIdMapT = std::unordered_map<std::string, int>;

  // Reads the next record from the specified input file. Returns the number of fields read, or 0
  // if the end of the file is reached.
  int nextRecord(std::ifstream& file) {
    assert(file);
    std::string lineStr;
    getline(file, lineStr);
    std::istringstream line(lineStr);
    if (!file)
      return 0;

    // Read the line token by token.
    int i = 0;
    while (!getline(line, currentRecord[i++], ',').eof()) {}
    return i;
  }

  std::ifstream verticesFile; // The input file containing the vertex records.
  std::ifstream edgesFile;    // The input file containing the edge records.

  RecordT currentRecord;  // A string array storing the current vertex/edge record.
  EdgeRecord currentEdge; // The current edge record.
  bool returnReverse;     // Indicates that next edge to be returned is the reverse of the current.

  std::vector<LatLng> latLngs; // The LatLng values of the vertices.
  VertexIdMapT origToNewIds;   // Maps original vertex IDs to new sequential IDs.
  int nextVertexId;            // The next free vertex ID.
};

// Returns the value of the LatLng attribute for the current vertex.
template <>
inline LatLngAttribute::Type XatfImporter::getValue<LatLngAttribute>() const {
  assert(!latLngs.empty());
  return latLngs.back();
}

// Returns the value of the travel time attribute for the current edge.
template <>
inline TravelTimeAttribute::Type XatfImporter::getValue<TravelTimeAttribute>() const {
  assert(currentEdge.xatfRoadCategory != XatfRoadCategory::UNUSED);
  const int cat = static_cast<int>(currentEdge.xatfRoadCategory);
  return currentEdge.xatfRoadCategory != XatfRoadCategory::FERRY
      ? std::round(36.0 * currentEdge.length / SPEED_LIMITS[cat])
      : currentEdge.length * TRAVEL_TIME_FACTOR;
}

// Returns the value of the length attribute for the current edge.
template <>
inline LengthAttribute::Type XatfImporter::getValue<LengthAttribute>() const {
  assert(currentEdge.xatfRoadCategory != XatfRoadCategory::UNUSED);
  return currentEdge.xatfRoadCategory != XatfRoadCategory::FERRY ? currentEdge.length
      : std::round(latLngs[edgeTail()].getGreatCircleDistanceTo(latLngs[edgeHead()]));
}

// Returns the value of the XATF road category attribute for the current edge.
template <>
inline XatfRoadCategoryAttribute::Type XatfImporter::getValue<XatfRoadCategoryAttribute>() const {
  return currentEdge.xatfRoadCategory;
}

// Returns the value of the capacity attribute for the current edge.
template <>
inline CapacityAttribute::Type XatfImporter::getValue<CapacityAttribute>() const {
  // Due to the lack of data, the capacity is computed from the number of lanes and speed limit of
  // the current edge, using the well-known two-second rule to determine a safe following distance.
  // Hence, capacity = numLanes * speedLimit / (2s * speedLimit + vehicleLength). The strange
  // factors in the expression below are due to unit conversions and avoiding divisions.
  assert(currentEdge.xatfRoadCategory != XatfRoadCategory::UNUSED);
  const int cat = static_cast<int>(currentEdge.xatfRoadCategory);
  return currentEdge.xatfRoadCategory == XatfRoadCategory::FERRY ? 1800 : std::round(
      9000.0 * NUM_LANES[cat] * SPEED_LIMITS[cat] / (5 * SPEED_LIMITS[cat] + 9 * VEHICLE_LENGTH));
}
