#pragma once

#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <unordered_map>

#include <rapidxml.hpp>

#include "DataStructures/Geometry/CoordinateTransformation.h"
#include "DataStructures/Geometry/LatLng.h"
#include "DataStructures/Geometry/Point.h"
#include "DataStructures/Graph/Attributes/CapacityAttribute.h"
#include "DataStructures/Graph/Attributes/CoordinateAttribute.h"
#include "DataStructures/Graph/Attributes/EdgeIdAttribute.h"
#include "DataStructures/Graph/Attributes/FreeFlowSpeedAttribute.h"
#include "DataStructures/Graph/Attributes/LatLngAttribute.h"
#include "DataStructures/Graph/Attributes/LengthAttribute.h"
#include "DataStructures/Graph/Attributes/NumLanesAttribute.h"
#include "DataStructures/Graph/Attributes/TravelTimeAttribute.h"
#include "Tools/Constants.h"
#include "Tools/DateHelpers.h"
#include "Tools/LexicalCast.h"
#include "Tools/StringHelpers.h"

// An importer to read graphs in MATSim network file format. The transport system (e.g., car, ride,
// freight, pt) whose network is to be imported must be passed to the constructor. First, the Graph
// class repeatedly calls nextVertex to read the next vertex from disk and fetches various vertex
// attributes. Then, it repeatedly calls nextEdge to read the next edge from disk and fetches
// various edge attributes.
class MatSimImporter {
 public:
  MatSimImporter(const std::string& system, const int epsgCode)
      : transportSystem(system),
        trans(epsgCode, CoordinateTransformation::WGS_84),
        analysisPeriod(1) {
    currentVertex.id = -1;
  }

  // Copy constructor.
  MatSimImporter(const MatSimImporter& im)
    : transportSystem(im.transportSystem), trans(im.trans), analysisPeriod(1) {
    currentVertex.id = -1;
  }

  // Opens the input file(s) and reads the header line(s).
  void init(const std::string& filename) {
    const auto fileSize = std::filesystem::file_size(filename + ".xml");
    networkBuffer.reset(new char[fileSize + 1]);
    std::ifstream networkFile(filename + ".xml");
    assert(networkFile.good());
    networkFile.read(networkBuffer.get(), fileSize);
    assert(networkFile.good());
    assert(networkFile.get() == std::ifstream::traits_type::eof());
    networkBuffer[fileSize] = '\0';
    networkDocument.parse<rapidxml::parse_default>(networkBuffer.get());

    const auto networkElem = networkDocument.first_node();
    assert(networkElem != nullptr);
    assert(stringEq(networkElem->name(), "network"));
    const auto nodesElem = networkElem->first_node();
    assert(nodesElem != nullptr);
    assert(stringEq(nodesElem->name(), "nodes"));
    const auto linksElem = nodesElem->next_sibling();
    assert(linksElem != nullptr);
    assert(stringEq(linksElem->name(), "links"));

    const auto capperiodAttr = linksElem->first_attribute();
    assert(capperiodAttr != nullptr);
    if (stringEq(capperiodAttr->name(), "capperiod")) {
      std::string capperiodVal(capperiodAttr->value());
      analysisPeriod = parseTime(capperiodVal) / 3600.0;
    }

    currentVertexElem = nodesElem->first_node();
    currentEdgeElem = linksElem->first_node();
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
    if (currentVertexElem == nullptr)
      return false;
    assert(stringEq(currentVertexElem->name(), "node"));

    const auto idAttr = currentVertexElem->first_attribute();
    assert(idAttr != nullptr);
    assert(stringEq(idAttr->name(), "id"));
    assert(origIdToSeqId.find(idAttr->value()) == origIdToSeqId.end());
    origIdToSeqId[idAttr->value()] = ++currentVertex.id;

    const auto xAttr = idAttr->next_attribute();
    assert(xAttr != nullptr);
    assert(stringEq(xAttr->name(), "x"));
    currentVertex.coordinate.x() = std::round(lexicalCast<double>(xAttr->value()));

    const auto yAttr = xAttr->next_attribute();
    assert(yAttr != nullptr);
    assert(stringEq(yAttr->name(), "y"));
    currentVertex.coordinate.y() = std::round(lexicalCast<double>(yAttr->value()));

    double lng, lat;
    trans.forward(currentVertex.coordinate.x(), currentVertex.coordinate.y(), lng, lat);
    currentVertex.latLng = {lat, lng};

    currentVertexElem = currentVertexElem->next_sibling();
    return true;
  }

  // Returns the ID of the current vertex. Vertices must have sequential IDs from 0 to n − 1.
  int vertexId() const {
    return currentVertex.id;
  }

  // Reads the next edge from disk. Returns false if there are no more edges.
  bool nextEdge() {
    while (currentEdgeElem != nullptr) {
      assert(stringEq(currentEdgeElem->name(), "link"));
      const auto modesAttr = currentEdgeElem->last_attribute();
      assert(modesAttr != nullptr);
      if (!stringEq(modesAttr->name(), "modes")) {
        currentEdgeElem = currentEdgeElem->next_sibling();
        continue;
      }
      const auto permittedSystems = modesAttr->value();

      // Set open to true if the selected transport system is permitted on the current edge.
      auto open = true;
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
      if (!open) {
        currentEdgeElem = currentEdgeElem->next_sibling();
        continue;
      }

      const auto idAttr = currentEdgeElem->first_attribute();
      assert(idAttr != nullptr);
      assert(stringEq(idAttr->name(), "id"));
      currentEdge.id = lexicalCast<int>(idAttr->value());

      const auto fromAttr = idAttr->next_attribute();
      assert(fromAttr != nullptr);
      assert(stringEq(fromAttr->name(), "from"));
      assert(origIdToSeqId.find(fromAttr->value()) != origIdToSeqId.end());
      currentEdge.tail = origIdToSeqId.at(fromAttr->value());

      const auto toAttr = fromAttr->next_attribute();
      assert(toAttr != nullptr);
      assert(stringEq(toAttr->name(), "to"));
      assert(origIdToSeqId.find(toAttr->value()) != origIdToSeqId.end());
      currentEdge.head = origIdToSeqId.at(toAttr->value());

      const auto lengthAttr = toAttr->next_attribute();
      assert(lengthAttr != nullptr);
      assert(stringEq(lengthAttr->name(), "length"));
      currentEdge.length = lexicalCast<double>(lengthAttr->value());

      const auto freespeedAttr = lengthAttr->next_attribute();
      assert(freespeedAttr != nullptr);
      assert(stringEq(freespeedAttr->name(), "freespeed"));
      currentEdge.freeFlowSpeed = std::round(lexicalCast<double>(freespeedAttr->value()) * 3.6);

      const auto capacityAttr = freespeedAttr->next_attribute();
      assert(capacityAttr != nullptr);
      assert(stringEq(capacityAttr->name(), "capacity"));
      currentEdge.capacity =
          std::round(lexicalCast<double>(capacityAttr->value()) / analysisPeriod);

      const auto permlanesAttr = capacityAttr->next_attribute();
      assert(permlanesAttr != nullptr);
      assert(stringEq(permlanesAttr->name(), "permlanes"));
      currentEdge.numLanes = lexicalCast<double>(permlanesAttr->value());

      currentEdgeElem = currentEdgeElem->next_sibling();
      return true;
    }

    return false;
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
  void close() {
    networkDocument.clear();
    networkBuffer.reset();
  }

 private:
  // A vertex record in MATSim network file format.
  struct VertexRecord {
    int id;
    Point coordinate;
    LatLng latLng;
  };

  // An edge record in MATSim network file format.
  struct EdgeRecord {
    int id;
    int tail;
    int head;
    double length;
    int freeFlowSpeed;
    int capacity;
    double numLanes;
  };

  using IdMap = std::unordered_map<std::string, int>;

  const std::string& transportSystem; // The system (car, freight) whose network is to be imported.
  CoordinateTransformation trans;     // Transformation from the input coordinate system to WGS84.

  std::unique_ptr<char[]> networkBuffer;    // A buffer in which the network file content is stored.
  rapidxml::xml_document<> networkDocument; // The XML-based network document.

  const rapidxml::xml_node<>* currentVertexElem; // The current vertex element in the XML document.
  const rapidxml::xml_node<>* currentEdgeElem;   // The current edge element in the XML document.
  VertexRecord currentVertex;                    // The current vertex record.
  EdgeRecord currentEdge;                        // The current edge record.

  IdMap origIdToSeqId;   // A map from original vertex IDs to sequential IDs.
  double analysisPeriod; // The analysis period in hours (capacity is given in vehicles/AP).
};

// Returns the value of the coordinate attribute for the current vertex.
template <>
inline CoordinateAttribute::Type MatSimImporter::getValue<CoordinateAttribute>() const {
  return currentVertex.coordinate;
}

// Returns the value of the edge ID attribute for the current edge.
template <>
inline EdgeIdAttribute::Type MatSimImporter::getValue<EdgeIdAttribute>() const {
  return currentEdge.id;
}

// Returns the value of the LatLng attribute for the current vertex.
template <>
inline LatLngAttribute::Type MatSimImporter::getValue<LatLngAttribute>() const {
  return currentVertex.latLng;
}

// Returns the value of the capacity attribute for the current edge.
template <>
inline CapacityAttribute::Type MatSimImporter::getValue<CapacityAttribute>() const {
  return currentEdge.capacity;
}

// Returns the value of the free-flow speed attribute for the current edge.
template <>
inline FreeFlowSpeedAttribute::Type MatSimImporter::getValue<FreeFlowSpeedAttribute>() const {
  return currentEdge.freeFlowSpeed;
}

// Returns the value of the length attribute for the current edge.
template <>
inline LengthAttribute::Type MatSimImporter::getValue<LengthAttribute>() const {
  return std::round(currentEdge.length);
}

// Returns the value of the number of lanes attribute for the current edge.
template <>
inline NumLanesAttribute::Type MatSimImporter::getValue<NumLanesAttribute>() const {
  return currentEdge.numLanes;
}

// Returns the value of the travel time attribute for the current edge.
template <>
inline TravelTimeAttribute::Type MatSimImporter::getValue<TravelTimeAttribute>() const {
  if (currentEdge.freeFlowSpeed != 0)
    return static_cast<int>(3.6 * currentEdge.length / currentEdge.freeFlowSpeed) * 10 + 10;
  else
    return INFTY;
}
