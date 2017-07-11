#pragma once

#include <cassert>
#include <vector>

#include "DataStructures/Geometry/Point.h"
#include "DataStructures/Geometry/Polygon.h"

namespace visum {

// A surface in the Visum surface model. It can consist of a number of polygonal regions, some with
// holes. A surface stores a list of simple polygons. Counterclockwise polygons add to the surface,
// clockwise polygons subtract from it. The polygons are added and subtracted in the order in which
// they appear in the list.
class Surface {
 public:
  // Returns true if p is inside the boundary of this surface.
  bool contains(const Point& p) const {
    for (auto face = faces.rbegin(); face != faces.rend(); ++face)
      if (face->contains(p))
        return face->orientation() > 0;
    return false;
  }

  // Appends a face to this surface.
  void addFace(const Polygon& face) {
    assert(face.simple());
    faces.push_back(face);
  }

 private:
  std::vector<Polygon> faces; // The polygons composing this surface.
};

}
