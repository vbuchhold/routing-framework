#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <string>
#include <unordered_map>
#include <utility>

#include <csv.h>

#include "DataStructures/Geometry/CoordinateTransformation.h"
#include "DataStructures/Geometry/LatLng.h"
#include "DataStructures/Geometry/Point.h"
#include "DataStructures/Graph/Attributes/CapacityAttribute.h"
#include "DataStructures/Graph/Attributes/CoordinateAttribute.h"
#include "DataStructures/Graph/Attributes/FreeFlowSpeedAttribute.h"
#include "DataStructures/Graph/Attributes/LatLngAttribute.h"
#include "DataStructures/Graph/Attributes/LengthAttribute.h"
#include "DataStructures/Graph/Attributes/NumLanesAttribute.h"
#include "DataStructures/Graph/Attributes/TravelTimeAttribute.h"
#include "DataStructures/Graph/Attributes/VertexIdAttribute.h"
#include "Tools/Constants.h"
#include "Tools/LexicalCast.h"
#include "Tools/Math.h"
#include "Tools/StringHelpers.h"

// An importer to read graphs in Visum network file format. The transport system (e.g. pedestrian,
// car, bicycle) whose network is to be imported must be supplied in the constructor. First, the
// Graph class repeatedly calls nextVertex to read the next vertex from disk and fetches various
// vertex attributes. Then, it repeatedly calls nextEdge to read the next edge from disk and
// fetches various edge attributes.
class VisumImporter {
 public:
  // Constructs an importer to read the specified system's network.
  VisumImporter(const std::string& filename, const std::string& system,
                const int epsgCode, const double coordinatePrecision, const double analysisPeriod)
      : vertexReader(filename + "/KNOTEN.csv"),
        edgeReader(filename + "/STRECKE.csv"),
        interPointReader(filename + "/STRECKENPOLY.csv"),
        transportSystem(system),
        coordinatePrecision(coordinatePrecision),
        analysisPeriod(analysisPeriod),
        trans(epsgCode, CoordinateTransformation::WGS_84),
        nextVertexId(0) {
    assert(coordinatePrecision > 0);
    assert(analysisPeriod > 0);
  }

  // Copy constructor.
  VisumImporter(const VisumImporter& im)
      : vertexReader(im.vertexReader.get_truncated_file_name()),
        edgeReader(im.edgeReader.get_truncated_file_name()),
        interPointReader(im.interPointReader.get_truncated_file_name()),
        transportSystem(im.transportSystem),
        coordinatePrecision(im.coordinatePrecision),
        analysisPeriod(im.analysisPeriod),
        trans(im.trans),
        nextVertexId(0) {}

  // Opens the input file(s) and reads the header line(s).
  void init(const std::string& filename) {
    vertexReader.read_header(io::ignore_extra_column, "NR", "XKOORD", "YKOORD");
    edgeReader.read_header(
        io::ignore_extra_column, "VONKNOTNR", "NACHKNOTNR", "TYPNR", "VSYSSET",
        "LAENGE", "ANZFAHRSTREIFEN", "KAPIV", "V0IV");
    interPointReader.read_header(
        io::ignore_extra_column, "VONKNOTNR", "NACHKNOTNR", "INDEX", "XKOORD", "YKOORD");

    // Read the first normal record line from the intermediate point file.
    currentInterPoint.tail = INVALID_ID;
    int idx;
    interPointReader.read_row(
        currentInterPoint.tail, currentInterPoint.head, idx,
        currentInterPoint.easting, currentInterPoint.northing);
    assert(currentInterPoint.tail == INVALID_ID || idx == 1);

    // Read the maximum speed of the selected system for each of the 100 edge types.
    VisumFileReader<2> edgeTypeReader(filename + "/STRECKENTYP.csv");
    std::string upperCaseTS = transportSystem;
    toUpperCase(upperCaseTS);
    edgeTypeReader.read_header(io::ignore_extra_column, "NR", "VMAX-IVSYS(" + upperCaseTS + ")");
    auto i = 0;
    int id;
    char* maxSpeedField;
    for (; edgeTypeReader.read_row(id, maxSpeedField); ++i) {
      assert(id == i);
      assert(endsWith(maxSpeedField, "km/h"));
      substr(maxSpeedField, 0, std::strlen(maxSpeedField) - 4);
      maxSpeed[i] = lexicalCast<int>(maxSpeedField);
      assert(maxSpeed[i] >= 0);
    }
    assert(i == 100);
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
    double easting, northing;
    if (!vertexReader.read_row(currentVertex.id, easting, northing))
      return false;

    assert(origToNewId.find(currentVertex.id) == origToNewId.end());
    origToNewId[currentVertex.id] = nextVertexId++;

    currentVertex.coordinate.x() = std::round(easting * coordinatePrecision);
    currentVertex.coordinate.y() = std::round(northing * coordinatePrecision);

    double lng, lat;
    trans.forward(
        currentVertex.coordinate.x() / coordinatePrecision,
        currentVertex.coordinate.y() / coordinatePrecision, lng, lat);
    currentVertex.latLng = {toDegrees(lat), toDegrees(lng)};
    return true;
  }

  // Returns the ID of the current vertex. Vertices must have sequential IDs from 0 to n − 1.
  int vertexId() const {
    return nextVertexId - 1;
  }

  // Reads the next edge from disk. Returns false if there are no more edges.
  bool nextEdge() {
    int edgeType;
    char *permittedSystems = nullptr, *lengthField, *speedField;
    bool open;

    do {
      if (!edgeReader.read_row(currentEdge.tail, currentEdge.head, edgeType, permittedSystems,
          lengthField, currentEdge.numLanes, currentEdge.capacity, speedField)) {
        assert(currentInterPoint.tail == INVALID_ID);
        return false;
      }

      // Set open to true if the selected transport system is permitted on the current edge.
      open = true;
      auto j = 0;
      for (auto i = 0; permittedSystems[i] != '\0'; ++i)
        if (permittedSystems[i] == ',') {
          if (open && transportSystem[j] == '\0')
            break;
          open = true;
          j = 0;
        } else if (open) {
          open = permittedSystems[i] == transportSystem[j++];
        }
      open &= transportSystem[j] == '\0';

      assert(edgeType >= 0); assert(edgeType < 100);
      assert(currentEdge.numLanes >= 0);
      assert(currentEdge.capacity >= 0);

      assert(endsWith(speedField, "km/h"));
      substr(speedField, 0, std::strlen(speedField) - 4);
      currentEdge.freeFlowSpeed = std::min(lexicalCast<int>(speedField), maxSpeed[edgeType]);
      assert(currentEdge.freeFlowSpeed >= 0);

      // Read the polyline of the current edge (if it has one).
      if ((currentEdge.tail == currentInterPoint.tail &&
          currentEdge.head == currentInterPoint.head) ||
          (currentEdge.tail == currentInterPoint.head &&
          currentEdge.head == currentInterPoint.tail)) {
        assert(origToNewId.find(currentInterPoint.tail) != origToNewId.end());
        assert(origToNewId.find(currentInterPoint.head) != origToNewId.end());
        currentPolyline.clear();
        currentPolylineEndpoints.first = origToNewId.at(currentInterPoint.tail);
        currentPolylineEndpoints.second = origToNewId.at(currentInterPoint.head);
        auto idx = 1;
        while ((currentEdge.tail == currentInterPoint.tail &&
               currentEdge.head == currentInterPoint.head) ||
               (currentEdge.tail == currentInterPoint.head &&
               currentEdge.head == currentInterPoint.tail)) {
          assert(currentInterPoint.tail >= 0);
          assert(currentInterPoint.head >= 0);
          assert(idx == currentPolyline.size() + 1);
          double lng, lat;
          trans.forward(currentInterPoint.easting, currentInterPoint.northing, lng, lat);
          currentPolyline.emplace_back(toDegrees(lat), toDegrees(lng));
          currentInterPoint.tail = INVALID_ID;
          interPointReader.read_row(
              currentInterPoint.tail, currentInterPoint.head, idx,
              currentInterPoint.easting, currentInterPoint.northing);
        }
        assert(currentInterPoint.tail == INVALID_ID || idx == 1);
      }
    } while (!open || currentEdge.capacity == 0 || currentEdge.freeFlowSpeed == 0);

    assert(origToNewId.find(currentEdge.tail) != origToNewId.end());
    assert(origToNewId.find(currentEdge.head) != origToNewId.end());
    currentEdge.tail = origToNewId.at(currentEdge.tail);
    currentEdge.head = origToNewId.at(currentEdge.head);

    assert(endsWith(lengthField, "km"));
    substr(lengthField, 0, std::strlen(lengthField) - 2);
    currentEdge.length = std::round(lexicalCast<double>(lengthField) * 1000);
    assert(currentEdge.length >= 0);
    return true;
  }

  // Returns the tail vertex of the current edge.
  int edgeTail() const {
    return currentEdge.tail;
  }

  // Returns the head vertex of the current edge.
  int edgeHead() const {
    return currentEdge.head;
  }

  // Returns the value of the specified attribute for the current vertex/edge, or the attribute's
  // default value if the attribute is not part of the file format.
  template <typename Attr>
  typename Attr::Type getValue() const {
    return Attr::defaultValue();
  }

  // Closes the input file(s).
  void close() { /* do nothing */ }

 private:
  // The CSV dialect used by Visum.
  template <int numFields>
  using VisumFileReader = io::CSVReader<numFields, io::trim_chars<>, io::no_quote_escape<';'>>;

  // A vertex record in Visum network file format.
  struct VertexRecord {
    int id;
    Point coordinate;
    LatLng latLng;
  };

  // An edge record in Visum network file format.
  struct EdgeRecord {
    int tail;
    int head;
    int length;
    int numLanes;
    int capacity;
    int freeFlowSpeed;
  };

  // An intermediate point record in the Visum network file format.
  struct InterPointRecord {
    int tail;
    int head;
    double easting;
    double northing;
  };

  using IdMap = std::unordered_map<int, int>;

  VisumFileReader<3> vertexReader;     // The CSV file that contains the vertex records.
  VisumFileReader<8> edgeReader;       // The CSV file that contains the edge records.
  VisumFileReader<5> interPointReader; // The CSV file that contains the intermediate point records.
  const std::string transportSystem;   // The system (car, bicycle) whose network is to be imported.
  const double coordinatePrecision;    // The number of digits to consider after the decimal point.
  const double analysisPeriod;            // The analysis period in hours (capacity is in vehicles/AP).
  CoordinateTransformation trans;      // Transformation from the input coordinate system to WGS84.

  int maxSpeed[100]; // The maximum speed of the selected system for each of the 100 edge types.
  IdMap origToNewId; // A map from original vertex IDs to new sequential IDs.
  int nextVertexId;  // The next free vertex ID.

  VertexRecord currentVertex;                   // The vertex read by the last call of nextVertex.
  EdgeRecord currentEdge;                       // The edge read by the last call of nextEdge.
  InterPointRecord currentInterPoint;           // The last intermediate point record read.
  std::vector<LatLng> currentPolyline;          // The last edge polyline read.
  std::pair<int, int> currentPolylineEndpoints; // The endpoints of the last edge polyline read.
};

// Returns the value of the coordinate attribute for the current vertex.
template <>
inline CoordinateAttribute::Type VisumImporter::getValue<CoordinateAttribute>() const {
  return currentVertex.coordinate;
}

// Returns the value of the LatLng attribute for the current vertex.
template <>
inline LatLngAttribute::Type VisumImporter::getValue<LatLngAttribute>() const {
  return currentVertex.latLng;
}

// Returns the value of the capacity attribute for the current edge.
template <>
inline CapacityAttribute::Type VisumImporter::getValue<CapacityAttribute>() const {
  return std::round(currentEdge.capacity / analysisPeriod);
}

// Returns the value of the free-flow speed attribute for the current edge.
template <>
inline FreeFlowSpeedAttribute::Type VisumImporter::getValue<FreeFlowSpeedAttribute>() const {
  return currentEdge.freeFlowSpeed;
}

// Returns the value of the length attribute for the current edge.
template <>
inline LengthAttribute::Type VisumImporter::getValue<LengthAttribute>() const {
  return currentEdge.length;
}

// Returns the value of the number of lanes attribute for the current edge.
template <>
inline NumLanesAttribute::Type VisumImporter::getValue<NumLanesAttribute>() const {
  return currentEdge.numLanes;
}

// Returns the value of the road geometry attribute for the current edge.
template <>
inline RoadGeometryAttribute::Type VisumImporter::getValue<RoadGeometryAttribute>() const {
  if (currentEdge.tail == currentPolylineEndpoints.first &&
      currentEdge.head == currentPolylineEndpoints.second)
    return currentPolyline;
  else if (currentEdge.tail == currentPolylineEndpoints.second &&
           currentEdge.head == currentPolylineEndpoints.first)
    return {currentPolyline.rbegin(), currentPolyline.rend()};
  return RoadGeometryAttribute::defaultValue();
}

// Returns the value of the travel time attribute for the current edge.
template <>
inline TravelTimeAttribute::Type VisumImporter::getValue<TravelTimeAttribute>() const {
  return std::round(36.0 * currentEdge.length / currentEdge.freeFlowSpeed);
}

// Returns the value of the vertex ID attribute for the current vertex.
template <>
inline VertexIdAttribute::Type VisumImporter::getValue<VertexIdAttribute>() const {
  return currentVertex.id;
}
