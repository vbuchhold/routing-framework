#pragma once

#include <cassert>
#include <cstring>
#include <fstream>
#include <ios>
#include <iterator>
#include <string>
#include <vector>

#include "DataStructures/Geometry/LatLng.h"
#include "DataStructures/Geometry/Point.h"
#include "DataStructures/Geometry/Polygon.h"
#include "DataStructures/Geometry/Rectangle.h"
#include "Tools/LexicalCast.h"
#include "Tools/StringHelpers.h"

// An area can consist of a number of polygonal regions, some with holes. It stores a list of simple
// polygons. Counterclockwise polygons add to the area, clockwise polygons subtract from it. The
// polygons are added and subtracted in the order in which they appear in the list. Areas can be
// read from and written to OSM POLY files.
class Area {
 public:
  // Iterators referring to faces of the area.
  using Iterator = std::vector<Polygon>::iterator;
  using ConstIterator = std::vector<Polygon>::const_iterator;

  // Returns an iterator referring to the first face of the area.
  ConstIterator begin() const {
    return faces.begin();
  }

  // Returns an iterator which is the past-the-end value for the area.
  ConstIterator end() const {
    return faces.end();
  }

  // Returns true if p is inside the boundary of this area.
  bool contains(const Point& p) const {
    for (auto face = faces.rbegin(); face != faces.rend(); ++face)
      if (face->contains(p))
        return face->orientation() > 0;
    return false;
  }

  // Returns the bounding box containing all vertices of this area.
  Rectangle boundingBox() const {
    int firstPositiveFace = 0;
    while (firstPositiveFace < faces.size() && faces[firstPositiveFace].orientation() < 0)
      ++firstPositiveFace;
    if (firstPositiveFace == faces.size())
      return Rectangle();
    Rectangle box(faces[firstPositiveFace].begin(), faces[firstPositiveFace].end());
    for (int i = 1; i < faces.size(); ++i)
      if (faces[i].orientation() > 0)
        box.extend(faces[i].begin(), faces[i].end());
    return box;
  }

  // Combines a polygon with this area. If the polygon is counterclockwise, it adds to the area.
  // Otherwise, it subtracts from the area.
  void combine(const Polygon& face) {
    assert(face.simple());
    faces.push_back(face);
  }

  // Adds a polygon to this area.
  void add(const Polygon& face) {
    assert(face.simple());
    if (face.orientation() > 0)
      faces.emplace_back(face.begin(), face.end());
    else {
      faces.emplace_back(face.rbegin(), face.rend());
    }
  }

  // Subtracts a polygon from this area.
  void subtract(const Polygon& face) {
    assert(face.simple());
    if (face.orientation() < 0)
      faces.emplace_back(face.begin(), face.end());
    else
      faces.emplace_back(face.rbegin(), face.rend());
  }

  // Reads an area from an OSM POLY file.
  void importFromOsmPolyFile(const std::string& filename) {
    faces.clear();
    std::ifstream in(filename);
    assert(in.good());
    std::string line;
    getline(in, line);
    bool inPolygon = false;
    bool addToArea = false;
    while (getline(in, line)) {
      trim(line);
      assert(!line.empty());
      if (!inPolygon && line == "END") {
        // End-of-file occurred.
        return;
      } else if (!inPolygon) {
        // Section header occurred.
        inPolygon = true;
        addToArea = line[0] != '!';
        faces.emplace_back();
      } else if (inPolygon && line == "END") {
        // End-of-section occurred.
        inPolygon = false;
        assert(faces.back().size() >= 2);
        if (faces.back()[1] == faces.back().back()) {
          faces.back().removeBack();
          faces.back().remove(faces.back().begin());
        } else if (faces.back()[0] == faces.back().back()) {
          faces.back().removeBack();
        }
        if ((faces.back().orientation() > 0) != addToArea)
          faces.back().reverseOrientation();
        assert(faces.back().simple());
      } else if (inPolygon) {
        // Vertex record occurred.
        const char* const lng = std::strtok(&line[0], " \f\n\r\t\v");
        const char* const lat = std::strtok(nullptr, " \f\n\r\t\v");
        assert(lng != nullptr);
        assert(lat != nullptr);
        assert(std::strtok(nullptr, " \f\n\r\t\v") == nullptr);
        LatLng latLng(lexicalCast<double>(lat), lexicalCast<double>(lng));
        faces.back().add({latLng.longitude(), latLng.latitude()});
      }
    }
    assert(false);
  }

  // Writes an area to an OSM POLY file.
  void exportToOsmPolyFile(const std::string& filename) const {
    std::ofstream out(filename);
    assert(out.good());
    out << std::uppercase << std::scientific;
    out << "none\n";
    for (int i = 0; i < faces.size(); ++i) {
      if (faces[i].orientation() < 0)
        out << '!';
      out << i + 1 << '\n';
      for (const auto& v : faces[i]) {
        LatLng latLng(v.getY(), v.getX());
        out << "  " << latLng.lngInDeg() << "  " << latLng.latInDeg() << '\n';
      }
      if (faces[i].size() > 0) {
        LatLng latLng(faces[i][0].getY(), faces[i][0].getX());
        out << "  " << latLng.lngInDeg() << "  " << latLng.latInDeg() << '\n';
      }
      out << "END\n";
    }
    out << "END\n";
  }

 private:
  std::vector<Polygon> faces; // The polygons composing this area.
};
