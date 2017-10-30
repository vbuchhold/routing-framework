#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "DataStructures/Geometry/LatLng.h"
#include "DataStructures/Graph/Attributes/CapacityAttribute.h"
#include "DataStructures/Graph/Attributes/LatLngAttribute.h"
#include "DataStructures/Graph/Attributes/LengthAttribute.h"
#include "DataStructures/Graph/Attributes/NumLanesAttribute.h"
#include "DataStructures/Graph/Attributes/TravelTimeAttribute.h"
#include "DataStructures/Graph/Attributes/XatfRoadCategoryAttribute.h"
#include "Tools/LexicalCast.h"

// An importer for reading graphs in XATF file format. First, the Graph class repeatedly calls
// nextVertex to read the next vertex from disk and fetches various vertex attributes. Then it
// repeatedly calls nextEdge to read the next edge from disk and fetches various edge attributes.
class XatfImporter {
 public:
  // Constructs an importer to read graphs in XATF file format.
  XatfImporter() : returnReverse(false), nextVertexId(0) {}

  // Opens the input file(s) and reads the header line(s).
  void init(const std::string& filename) {
    vertexFile.open(filename + ".nfx");
    edgeFile.open(filename + ".sfx");
    assert(vertexFile.good());
    assert(edgeFile.good());
  }

  // Returns the number of vertices in the graph, or 0 if the number is not yet known.
  int numVertices() const {
    return 0;
  }

  // Returns the number of edges in the graph, or 0 if the number is not yet known.
  int numEdges() const {
    return 0;
  }

  // Reads the next vertex from disk. Returns false if there are no more vertices.
  bool nextVertex() {
    int fieldsRead = nextRecord(vertexFile);
    if (fieldsRead == 0)
      return false;
    assert(fieldsRead == NUM_VERTEX_FIELDS);

    // Map the original vertex ID to a sequential ID.
    currentVertex.id = lexicalCast<int>(currentRecord[VERTEX_ID]);
    assert(origToNewIds.find(currentVertex.id) == origToNewIds.end());
    origToNewIds[currentVertex.id] = nextVertexId;
    currentVertex.id = nextVertexId++;

    // Store the LatLng of the vertex (possibly needed to compute edge lengths).
    const int lat = LatLng::PRECISION / 100000 * lexicalCast<int>(currentRecord[LATITUDE]);
    const int lng = LatLng::PRECISION / 100000 * lexicalCast<int>(currentRecord[LONGITUDE]);
    latLngs.emplace_back(lat, lng);
    assert(latLngs.size() == nextVertexId);
    return true;
  }

  // Returns the ID of the current vertex. Vertices must have sequential IDs from 0 to n − 1.
  int vertexId() const {
    return currentVertex.id;
  }

  // Reads the next edge from disk. Returns false if there are no more edges.
  bool nextEdge() {
    if (returnReverse) {
      std::swap(currentEdge.edgeTail, currentEdge.edgeHead);
      returnReverse = false;
      return true;
    }

    int fieldsRead = nextRecord(edgeFile);
    if (fieldsRead == 0)
      return false;
    assert(fieldsRead == NUM_EDGE_FIELDS);

    currentEdge.edgeTail = lexicalCast<int>(currentRecord[EDGE_TAIL]);
    currentEdge.edgeHead = lexicalCast<int>(currentRecord[EDGE_HEAD]);
    assert(origToNewIds.find(currentEdge.edgeTail) != origToNewIds.end());
    assert(origToNewIds.find(currentEdge.edgeHead) != origToNewIds.end());
    currentEdge.edgeTail = origToNewIds.at(currentEdge.edgeTail);
    currentEdge.edgeHead = origToNewIds.at(currentEdge.edgeHead);

    currentEdge.length = lexicalCast<int>(currentRecord[LENGTH]);
    currentEdge.xatfRoadCategory = lexicalCast<XatfRoadCategory>(currentRecord[XATF_ROAD_CATEGORY]);
    assert(currentEdge.length >= 0);

    // Handle edge direction. We regard closed edges as bidirected to obtain a larger graph.
    Direction dir = lexicalCast<Direction>(currentRecord[DIRECTION]);
    assert(dir == BIDIRECTED || dir == FORWARD || dir == REVERSE || dir == CLOSED);
    returnReverse = dir == BIDIRECTED || dir == CLOSED;
    if (dir == REVERSE)
      std::swap(currentEdge.edgeTail, currentEdge.edgeHead);
    return true;
  }

  // Returns the tail vertex of the current edge.
  int edgeTail() const {
    return currentEdge.edgeTail;
  }

  // Returns the head vertex of the current edge.
  int edgeHead() const {
    return currentEdge.edgeHead;
  }

  // Returns the value of the specified attribute for the current vertex/edge, or the attribute's
  // default value if the attribute is not part of the file format.
  template <typename Attr>
  typename Attr::Type getValue() const {
    return Attr::defaultValue();
  }

  // Closes the input file(s).
  void close() {
    vertexFile.close();
    edgeFile.close();
  }

 private:
  // A vertex record in XATF file format.
  struct VertexRecord {
    int id;
  };

  // An edge record in XATF file format.
  struct EdgeRecord {
    int edgeTail;
    int edgeHead;
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

  // This enum specifies possible values of the "direction" edge field.
  enum Direction {
    BIDIRECTED = 0,
    FORWARD = 1,
    REVERSE = 2,
    CLOSED  = 3,
  };

  // The number of fields in each vertex/edge record.
  static constexpr int NUM_VERTEX_FIELDS = 5;
  static constexpr int NUM_EDGE_FIELDS   = 10;
  static constexpr int MAX_NUM_FIELDS = std::max(NUM_VERTEX_FIELDS, NUM_EDGE_FIELDS);

  // The free-flow speed in km/h for each XATF road category.
  static constexpr const int (&FREE_FLOW_SPEED)[] = {
    -1,
    130, // Motorway, fast
    120, // Motorway, medium
    110, // Motorway, slow
    100, // National road, fast
    90,  // National road, medium
    80,  // National road, slow
    70,  // Regional road, fast
    60,  // Regional road, medium
    50,  // Regional road, slow
    40,  // Urban street, fast
    30,  // Urban street, medium
    20,  // Urban street, slow
    -1,  // Ferry
    -1,  // Unused
    10,  // Forest road
  };

  // The number of lanes for each XATF road category.
  static constexpr const int (&NUM_LANES)[] = {
    -1,
    2,  // Motorway, fast
    2,  // Motorway, medium
    2,  // Motorway, slow
    2,  // National road, fast
    1,  // National road, medium
    1,  // National road, slow
    1,  // Regional road, fast
    1,  // Regional road, medium
    1,  // Regional road, slow
    1,  // Urban street, fast
    1,  // Urban street, medium
    1,  // Urban street, slow
    -1, // Ferry
    -1, // Unused
    1,  // Forest road
  };

  static constexpr int VEHICLE_LENGTH = 5; // The average vehicle length in meters.

  // Reads the next record from the specified input file. Returns the number of fields read, or 0
  // if the end of the file has been reached.
  int nextRecord(std::ifstream& file) {
    std::string lineStr;
    getline(file, lineStr);
    std::istringstream line(lineStr);
    if (!file)
      return 0;

    // Read the line field by field.
    int i = 0;
    while (!getline(line, currentRecord[i++], ',').eof()) {}
    return i;
  }

  using Record = std::array<std::string, MAX_NUM_FIELDS>;
  using IdMap = std::unordered_map<int, int>;

  std::ifstream vertexFile; // The input file containing the vertex records.
  std::ifstream edgeFile;   // The input file containing the edge records.

  Record currentRecord;       // A string array storing the current vertex/edge record.
  VertexRecord currentVertex; // The vertex record read by the last call of nextVertex.
  EdgeRecord currentEdge;     // The edge record read by the last call of nextEdge.
  bool returnReverse;         // Indicates that next edge to be returned is reverse of the current.

  std::vector<LatLng> latLngs; // The LatLngs of the vertices.
  IdMap origToNewIds;          // A map from original vertex IDs to new sequential IDs.
  int nextVertexId;            // The next free vertex ID.
};

// Returns the value of the LatLng attribute for the current vertex.
template <>
inline LatLngAttribute::Type XatfImporter::getValue<LatLngAttribute>() const {
  assert(!latLngs.empty());
  return latLngs.back();
}

// Returns the value of the length attribute for the current edge.
template <>
inline LengthAttribute::Type XatfImporter::getValue<LengthAttribute>() const {
  assert(currentEdge.xatfRoadCategory != XatfRoadCategory::UNUSED);
  if (currentEdge.xatfRoadCategory != XatfRoadCategory::FERRY)
    return currentEdge.length;
  else
    return std::round(latLngs[edgeTail()].getGreatCircleDistanceTo(latLngs[edgeHead()]));
}

// Returns the value of the XATF road category attribute for the current edge.
template <>
inline XatfRoadCategoryAttribute::Type XatfImporter::getValue<XatfRoadCategoryAttribute>() const {
  return currentEdge.xatfRoadCategory;
}

// Returns the value of the travel time attribute for the current edge.
template <>
inline TravelTimeAttribute::Type XatfImporter::getValue<TravelTimeAttribute>() const {
  assert(currentEdge.xatfRoadCategory != XatfRoadCategory::UNUSED);
  const int cat = static_cast<int>(currentEdge.xatfRoadCategory);
  if (currentEdge.xatfRoadCategory != XatfRoadCategory::FERRY)
    return std::round(36.0 * currentEdge.length / FREE_FLOW_SPEED[cat]);
  else
    return 600 * currentEdge.length;
}

// Returns the value of the number of lanes attribute for the current edge.
template <>
inline NumLanesAttribute::Type XatfImporter::getValue<NumLanesAttribute>() const {
  assert(currentEdge.xatfRoadCategory != XatfRoadCategory::UNUSED);
  return NUM_LANES[static_cast<int>(currentEdge.xatfRoadCategory)];
}

// Returns the value of the capacity attribute for the current edge.
template <>
inline CapacityAttribute::Type XatfImporter::getValue<CapacityAttribute>() const {
  // Due to lack of data, the capacity is computed from the number of lanes and free-flow speed of
  // the current edge, using the well-known two-second rule to determine a safe following distance.
  // Hence, capacity = numLanes * freeFlowSpeed / (2s * freeFlowSpeed + vehicleLength). The strange
  // factors in the expression below are due to unit conversions and avoiding divisions.
  assert(currentEdge.xatfRoadCategory != XatfRoadCategory::UNUSED);
  const int cat = static_cast<int>(currentEdge.xatfRoadCategory);
  if (currentEdge.xatfRoadCategory == XatfRoadCategory::FERRY)
    return 0;
  else
    return std::round(9000.0 *
        NUM_LANES[cat] * FREE_FLOW_SPEED[cat] / (5 * FREE_FLOW_SPEED[cat] + 9 * VEHICLE_LENGTH));
}
