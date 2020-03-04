#pragma once

#include <cassert>
#include <cmath>
#include <cstring>
#include <fstream>
#include <ios>
#include <limits>
#include <string>

#include "DataStructures/Geometry/LatLng.h"
#include "DataStructures/Graph/Attributes/LatLngAttribute.h"
#include "DataStructures/Graph/Attributes/LengthAttribute.h"
#include "DataStructures/Graph/Attributes/TravelTimeAttribute.h"
#include "Tools/Constants.h"

// An importer for reading graphs in DIMACS file format. First, the Graph class repeatedly calls
// nextVertex to read the next vertex from disk and fetches various vertex attributes. Then it
// repeatedly calls nextEdge to read the next edge from disk and fetches various edge attributes.
class DimacsImporter {
 public:
  // Constructs an importer to read graphs in DIMACS file format.
  DimacsImporter(const int distPrecision, const int timePrecision, const int coordinatePrecision)
      : distPrecision(distPrecision),
        timePrecision(timePrecision),
        coordinatePrecision(coordinatePrecision),
        vertexCount(0),
        edgeCount(0) {
    assert(distPrecision > 0);
    assert(timePrecision > 0);
    assert(coordinatePrecision > 0);
  }

  // Returns the number of vertices in the graph, or 0 if the number is not yet known.
  int numVertices() const {
    return vertexCount;
  }

  // Returns the number of edges in the graph, or 0 if the number is not yet known.
  int numEdges() const {
    return edgeCount;
  }

  // Opens the input file(s) and reads the header line(s).
  void init(const std::string& filename) {
    distGrFile.open(filename + ".dist.gr");
    if (distGrFile.is_open()) {
      while (distGrFile.peek() == std::char_traits<char>::to_int_type('c'))
        distGrFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      char header[4];
      distGrFile.read(header, 4) >> vertexCount >> edgeCount;
      distGrFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      assert(!distGrFile.fail());
      assert(std::memcmp(header, "p sp", 4) == 0);
      assert(vertexCount >= 0);
      assert(edgeCount >= 0);
    } else {
      currentEdge.dist = LengthAttribute::defaultValue();
    }

    timeGrFile.open(filename + ".time.gr");
    if (timeGrFile.is_open()) {
      while (timeGrFile.peek() == std::char_traits<char>::to_int_type('c'))
        timeGrFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      char header[4]; int vertexCount, edgeCount;
      timeGrFile.read(header, 4) >> vertexCount >> edgeCount;
      timeGrFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      assert(!timeGrFile.fail());
      assert(std::memcmp(header, "p sp", 4) == 0);
      assert(vertexCount >= 0);
      assert(edgeCount >= 0);

      if (distGrFile.is_open()) {
        assert(this->vertexCount == vertexCount);
        assert(this->edgeCount == edgeCount);
      } else {
        this->vertexCount = vertexCount;
        this->edgeCount = edgeCount;
      }
    } else {
      assert(distGrFile.is_open());
      currentEdge.time = TravelTimeAttribute::defaultValue();
    }

    coFile.open(filename + ".co");
    if (coFile.is_open()) {
      while (coFile.peek() == std::char_traits<char>::to_int_type('c'))
        coFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      char header[11]; int vertexCount;
      coFile.read(header, 11) >> vertexCount;
      coFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      assert(!coFile.fail());
      assert(std::memcmp(header, "p aux sp co", 11) == 0);
      assert(this->vertexCount == vertexCount);
    } else {
      currentVertex.id = 0;
      currentVertex.latLng = LatLngAttribute::defaultValue();
    }
  }

  // Reads the next vertex from disk. Returns false if there are no more vertices.
  bool nextVertex() {
    if (coFile.is_open()) {
      while (coFile.peek() == std::char_traits<char>::to_int_type('c'))
        coFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      char type; int lng, lat;
      coFile >> type >> currentVertex.id >> lng >> lat;
      coFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      if (!coFile.eof()) {
        assert(!coFile.fail());
        assert(type == 'v');
        assert(currentVertex.id > 0); assert(currentVertex.id <= vertexCount);
      }
      const auto conversionFactor = static_cast<double>(LatLng::PRECISION) / coordinatePrecision;
      currentVertex.latLng = {conversionFactor * lat, conversionFactor * lng};
      return !coFile.eof();
    } else {
      return ++currentVertex.id <= vertexCount;
    }
  }

  // Returns the ID of the current vertex. Vertices must have sequential IDs from 0 to n − 1.
  int vertexId() const {
    return currentVertex.id - 1;
  }

  // Reads the next edge from disk. Returns false if there are no more edges.
  bool nextEdge() {
    if (distGrFile.is_open()) {
      while (distGrFile.peek() == std::char_traits<char>::to_int_type('c'))
        distGrFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      char type;
      distGrFile >> type >> currentEdge.tail >> currentEdge.head >> currentEdge.dist;
      distGrFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      if (!distGrFile.eof()) {
        assert(!distGrFile.fail());
        assert(type == 'a');
      }
      currentEdge.dist = std::round(1.0 / distPrecision * currentEdge.dist);
      assert(currentEdge.dist >= 0); assert(currentEdge.dist < INFTY);
    }

    if (timeGrFile.is_open()) {
      while (timeGrFile.peek() == std::char_traits<char>::to_int_type('c'))
        timeGrFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      char type; int tail, head;
      timeGrFile >> type >> tail >> head >> currentEdge.time;
      timeGrFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      if (!timeGrFile.eof()) {
        assert(!timeGrFile.fail());
        assert(type == 'a');
      }
      currentEdge.time = std::round(10.0 / timePrecision * currentEdge.time);
      assert(currentEdge.time >= 0); assert(currentEdge.time < INFTY);

      if (distGrFile.is_open()) {
        assert(distGrFile.eof() == timeGrFile.eof());
        if (!distGrFile.eof()) {
          assert(currentEdge.tail == tail);
          assert(currentEdge.head == head);
        }
      } else {
        currentEdge.tail = tail;
        currentEdge.head = head;
      }
    }

    if (!(distGrFile.eof() || timeGrFile.eof())) {
      assert(currentEdge.tail > 0); assert(currentEdge.tail <= vertexCount);
      assert(currentEdge.head > 0); assert(currentEdge.head <= vertexCount);
      return true;
	}
    return false;
  }

  // Returns the tail vertex of the current edge.
  int edgeTail() const {
    return currentEdge.tail - 1;
  }

  // Returns the head vertex of the current edge.
  int edgeHead() const {
    return currentEdge.head - 1;
  }

  // Returns the value of the specified attribute for the current vertex/edge, or the attribute's
  // default value if the attribute is not part of the file format.
  template <typename Attr>
  typename Attr::Type getValue() const {
    return Attr::defaultValue();
  }

  // Closes the input file(s).
  void close() {
    distGrFile.close();
    timeGrFile.close();
    coFile.close();
  }

 private:
  // A vertex record in DIMACS file format.
  struct VertexRecord {
    int id;
    LatLng latLng;
  };

  // An edge record in DIMACS file format.
  struct EdgeRecord {
    int tail;
    int head;
    int dist;
    int time;
  };

  int distPrecision;       // Travel distances are given in 1/distPrecision meters.
  int timePrecision;       // Travel times are given in 1/timePrecision seconds.
  int coordinatePrecision; // Coordinates are given in 1/coordinatePrecision degrees.

  std::ifstream distGrFile; // The graph file that contains the travel distance metric.
  std::ifstream timeGrFile; // The graph file that contains the travel time metric.
  std::ifstream coFile;     // The auxiliary file that contains the coordinates.

  int vertexCount; // The number of vertices in the graph.
  int edgeCount;   // The number of edges in the graph.

  VertexRecord currentVertex; // The vertex record read by the last call of nextVertex.
  EdgeRecord currentEdge;     // The edge record read by the last call of nextEdge.
};

// Returns the value of the LatLng attribute for the current vertex.
template <>
inline LatLngAttribute::Type DimacsImporter::getValue<LatLngAttribute>() const {
  return currentVertex.latLng;
}

// Returns the value of the length attribute for the current edge.
template <>
inline LengthAttribute::Type DimacsImporter::getValue<LengthAttribute>() const {
  return currentEdge.dist;
}

// Returns the value of the travel time attribute for the current edge.
template <>
inline TravelTimeAttribute::Type DimacsImporter::getValue<TravelTimeAttribute>() const {
  return currentEdge.time;
}
