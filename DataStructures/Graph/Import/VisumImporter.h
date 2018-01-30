#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <string>
#include <unordered_map>

#include <csv.h>

#include "DataStructures/Geometry/CoordinateConversion.h"
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
#include "Tools/LexicalCast.h"
#include "Tools/StringHelpers.h"

// An importer to read graphs in Visum network file format. The transport system (e.g. pedestrian,
// bicycle, car) whose network is to be imported must be supplied in the constructor. First, the
// Graph class repeatedly calls nextVertex to read the next vertex from disk and fetches various
// vertex attributes. Then, it repeatedly calls nextEdge to read the next edge from disk and
// fetches various edge attributes.
class VisumImporter {
 public:
  // Constructs an importer to read the specified system's network.
  VisumImporter(const std::string& filename,
                const std::string& system, const int epsgCode, const int analysisPeriod)
      : vertexFile(filename + "/KNOTEN.csv"),
        edgeFile(filename + "/STRECKE.csv"),
        transportSystem(system),
        analysisPeriod(analysisPeriod),
        conversion(epsgCode),
        nextVertexId(0) {
    assert(analysisPeriod > 0);
  }

  // Copy constructor.
  VisumImporter(const VisumImporter& im)
      : vertexFile(im.vertexFile.get_truncated_file_name()),
        edgeFile(im.edgeFile.get_truncated_file_name()),
        transportSystem(im.transportSystem),
        analysisPeriod(im.analysisPeriod),
        conversion(im.conversion),
        nextVertexId(0) {}

  // Opens the input file(s) and reads the header line(s).
  void init(const std::string& filename) {
    vertexFile.read_header(io::ignore_extra_column, "NR", "XKOORD", "YKOORD");
    edgeFile.read_header(
        io::ignore_extra_column, "VONKNOTNR", "NACHKNOTNR", "TYPNR", "VSYSSET",
        "LAENGE", "ANZFAHRSTREIFEN", "KAPIV", "V0IV");

    // Read the maximum speed of the selected system for each of the 100 edge types.
    CsvDialect<2> edgeTypeFile(filename + "/STRECKENTYP.csv");
    std::string upperCaseTS = transportSystem;
    toUpperCase(upperCaseTS);
    edgeTypeFile.read_header(io::ignore_extra_column, "NR", "VMAX-IVSYS(" + upperCaseTS + ")");
    int i = 0;
    int id;
    char* maxSpeedField;
    for (; edgeTypeFile.read_row(id, maxSpeedField); ++i) {
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
    if (!vertexFile.read_row(currentVertex.id,
        currentVertex.coordinate.getX(), currentVertex.coordinate.getY()))
      return false;

    assert(origToNewIds.find(currentVertex.id) == origToNewIds.end());
    origToNewIds[currentVertex.id] = nextVertexId++;
    currentVertex.latLng = conversion.convert(currentVertex.coordinate);
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
      if (!edgeFile.read_row(currentEdge.edgeTail, currentEdge.edgeHead, edgeType,
          permittedSystems, lengthField, currentEdge.numLanes, currentEdge.capacity, speedField))
        return false;

      // Set open to true if the selected transport system is permitted on the current edge.
      open = true;
      int j = 0;
      for (int i = 0; permittedSystems[i] != '\0'; ++i)
        if (permittedSystems[i] == ',') {
          if (open && transportSystem[j] == '\0')
            break;
          open = true;
          j = 0;
        } else if (open) {
          open = permittedSystems[i] == transportSystem[j++];
        }
      open &= transportSystem[j] == '\0';
      if (!open)
        continue;

      assert(edgeType >= 0); assert(edgeType < 100);
      assert(currentEdge.numLanes >= 0);
      assert(currentEdge.capacity >= 0);

      assert(endsWith(speedField, "km/h"));
      substr(speedField, 0, std::strlen(speedField) - 4);
      currentEdge.freeFlowSpeed = std::min(lexicalCast<int>(speedField), maxSpeed[edgeType]);
      assert(currentEdge.freeFlowSpeed >= 0);
    } while (!open || currentEdge.capacity == 0 || currentEdge.freeFlowSpeed == 0);

    assert(origToNewIds.find(currentEdge.edgeTail) != origToNewIds.end());
    assert(origToNewIds.find(currentEdge.edgeHead) != origToNewIds.end());
    currentEdge.edgeTail = origToNewIds.at(currentEdge.edgeTail);
    currentEdge.edgeHead = origToNewIds.at(currentEdge.edgeHead);

    assert(endsWith(lengthField, "km"));
    substr(lengthField, 0, std::strlen(lengthField) - 2);
    currentEdge.length = std::round(lexicalCast<double>(lengthField) * 1000);
    assert(currentEdge.length >= 0);
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
  void close() { /* do nothing */ }

 private:
  // The CSV dialect used by Visum.
  template <int numFields>
  using CsvDialect = io::CSVReader<numFields, io::trim_chars<>, io::no_quote_escape<';'>>;

  // A vertex record in Visum network file format.
  struct VertexRecord {
    int id;
    Point coordinate;
    LatLng latLng;
  };

  // An edge record in Visum network file format.
  struct EdgeRecord {
    int edgeTail;
    int edgeHead;
    int length;
    int numLanes;
    int capacity;
    int freeFlowSpeed;
  };

  using IdMap = std::unordered_map<int, int>;

  CsvDialect<3> vertexFile;          // The CSV file containing the vertex records.
  CsvDialect<8> edgeFile;            // The CSV file containing the edge records.
  const std::string transportSystem; // The system (bicycle, car) whose network is to be imported.
  const int analysisPeriod;          // The analysis period in hours (capacity is in vehicles/AP).
  CoordinateConversion conversion;   // Conversion between the input coordinate system and WGS84.

  int maxSpeed[100];  // The maximum speed of the selected system for each of the 100 edge types.
  IdMap origToNewIds; // A map from original vertex IDs to new sequential IDs.
  int nextVertexId;   // The next free vertex ID.

  VertexRecord currentVertex; // The vertex record read by the last call of nextVertex.
  EdgeRecord currentEdge;     // The edge record read by the last call of nextEdge.
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
  return std::round(1.0 * currentEdge.capacity / analysisPeriod);
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
