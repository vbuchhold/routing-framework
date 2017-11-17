#pragma once

#include <cassert>
#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <csv.h>

#include "DataStructures/Geometry/Area.h"
#include "DataStructures/Geometry/Helpers.h"
#include "DataStructures/Geometry/Point.h"
#include "DataStructures/Geometry/Polygon.h"
#include "Tools/Constants.h"
#include "Tools/Workarounds.h"

namespace visum {

// The CSV dialect used by Visum.
template <int numFields>
using CsvDialect = io::CSVReader<numFields, io::trim_chars<>, io::no_quote_escape<';'>>;

// Follows the boundary of p until encountering a vertex a second time and splits p at the cycle
// found. The edges forming the cycle are stored in p2, the other edges in p1.
inline void splitPolygon(const Polygon& p, Polygon& p1, Polygon& p2) {
  int i, j = 0;
  bool cycleFound = false;
  for (i = 0; i < p.size(); ++i) {
    for (j = i + 1; j < p.size(); ++j)
      if (p[i] == p[j]) {
        cycleFound = true;
        break;
      }
    if (cycleFound)
      break;
  }
  p2 = {p.begin() + i, p.begin() + j};
  p1 = {p.begin() + j, p.end()};
  p1.add(p.begin(), p.begin() + i);
}

// Returns all zone polygons (and their IDs) from the specified Visum network file.
inline std::map<int, Area> readZonePolygonsFrom(const std::string& filename) {
  // Read all vertices from the Visum network file.
  std::unordered_map<int, Point> vertex;
  {
    int id, x, y;
    CsvDialect<3> vertexFile(filename + "/PUNKT.csv");
    vertexFile.read_header(io::ignore_no_column, "ID", "XKOORD", "YKOORD");
    while (vertexFile.read_row(id, x, y)) {
      assert(vertex.find(id) == vertex.end());
      vertex.emplace(id, Point(x, y));
    }
  }

  // Read all edges from the Visum network file.
  std::unordered_map<int, std::vector<Point>> edge;
  {
    int id, tail, head, edgeId = INVALID_ID, idx, x, y;
    CsvDialect<3> edgeFile(filename + "/KANTE.csv");
    CsvDialect<4> middleVertexFile(filename + "/ZWISCHENPUNKT.csv");
    edgeFile.read_header(io::ignore_no_column, "ID", "VONPUNKTID", "NACHPUNKTID");
    middleVertexFile.read_header(io::ignore_no_column, "KANTEID", "INDEX", "XKOORD", "YKOORD");
    middleVertexFile.read_row(edgeId, idx, x, y);
    while (edgeFile.read_row(id, tail, head)) {
      assert(id >= 0);
      assert(vertex.find(tail) != vertex.end());
      assert(vertex.find(head) != vertex.end());
      std::vector<Point>& currentEdge = edge[id];
      assert(currentEdge.empty());
      currentEdge.push_back(vertex[tail]);
      int prevIdx = 0;
      unused(prevIdx);
      while (edgeId == id) {
        edgeId = INVALID_ID;
        assert(idx == prevIdx + 1);
        const int len = currentEdge.size();
        Point p(x, y);
        if (len > 1 && p == currentEdge[len - 2])
          currentEdge.pop_back();
        else if (p != currentEdge[len - 1])
          currentEdge.push_back(p);
        prevIdx = idx;
        middleVertexFile.read_row(edgeId, idx, x, y);
      }
      const int len = currentEdge.size();
      Point p = vertex[head];
      if (len > 1 && p == currentEdge[len - 2])
        currentEdge.pop_back();
      else if (p != currentEdge[len - 1])
        currentEdge.push_back(p);
    }
    assert(edgeId == INVALID_ID);
    std::unordered_map<int, Point> tmp;
    vertex.swap(tmp);
  }

  // Read all polygons from the Visum network file.
  std::unordered_map<int, Polygon> polygon;
  {
    int id, idx, edgeId, dir;
    int prevId = INVALID_ID, prevIdx = INVALID_INDEX;
    unused(prevIdx);
    CsvDialect<4> polygonFile(filename + "/TEILFLAECHENELEMENT.csv");
    polygonFile.read_header(io::ignore_no_column, "TFLAECHEID", "INDEX", "KANTEID", "RICHTUNG");
    while (polygonFile.read_row(id, idx, edgeId, dir)) {
      assert(id >= 0);
      assert(edge.find(edgeId) != edge.end());
      assert(dir == 0 || dir == 1);
      Polygon& currentPolygon = polygon[id];
      const std::vector<Point>& currentEdge = edge[edgeId];
      if (id != prevId) {
        assert(currentPolygon.empty());
        assert(idx == 1);
        if (prevId != INVALID_ID) {
          const auto prevPolygon = polygon.find(prevId);
          assert(prevPolygon->second.front() == prevPolygon->second.back());
          prevPolygon->second.removeBack();
          if (!prevPolygon->second.simple()) {
            Polygon p1, p2;
            splitPolygon(prevPolygon->second, p1, p2);
            bool good = false;
            if (!p2.empty()) {
              good = true;
              for (int i = p1.size() - 1, j = 0; j < p1.size(); i = j++)
                for (int k = p2.size() - 1, l = 0; l < p2.size(); k = l++)
                  if (j > 1 || l > 1)
                    good &= !intersection(p1[i], p1[j], p2[k], p2[l]);
              const bool nested = p1.contains(p2.back()) || p2.contains(p1.back());
              good &= p1.simple() && p2.simple() &&
                  (p1.orientation() != p2.orientation()) == nested;
            }
            if (!good)
              polygon.erase(prevPolygon);
          }
        }
        if (dir == 0)
          currentPolygon.add(currentEdge.front());
        else
          currentPolygon.add(currentEdge.back());
      } else {
        assert(!currentPolygon.empty());
        assert(idx == prevIdx + 1);
      }
      if (dir == 0) {
        assert(currentPolygon.back() == currentEdge.front());
        currentPolygon.add(currentEdge.begin() + 1, currentEdge.end());
      } else {
        assert(currentPolygon.back() == currentEdge.back());
        currentPolygon.add(currentEdge.rbegin() + 1, currentEdge.rend());
      }
      prevId = id;
      prevIdx = idx;
    }
    assert(polygon[prevId].front() == polygon[prevId].back());
    polygon[prevId].removeBack();
    assert(polygon[prevId].simple());
    std::unordered_map<int, std::vector<Point>> tmp;
    edge.swap(tmp);
  }

  // Read all surfaces from the Visum network file.
  std::unordered_map<int, Area> surface;
  {
    int id, polygonId, hole;
    int prevId = INVALID_ID;
    bool broken = false;
    CsvDialect<3> surfaceFile(filename + "/FLAECHENELEMENT.csv");
    surfaceFile.read_header(io::ignore_no_column, "FLAECHEID", "TFLAECHEID", "ENKLAVE");
    while (surfaceFile.read_row(id, polygonId, hole)) {
      assert(id >= 0);
      assert(id >= prevId);
      assert(hole == 0 || hole == 1);
      assert(id == prevId || hole == 0);
      broken &= id == prevId;
      if (broken)
        continue;
      const auto currentPolygon = polygon.find(polygonId);
      if (currentPolygon == polygon.end()) {
        broken = true;
      } else if (currentPolygon->second.simple()) {
        if ((currentPolygon->second.orientation() > 0) == (hole == 0))
          surface[id].combine(currentPolygon->second);
        else
          broken = true;
      } else {
        Polygon p1, p2;
        splitPolygon(currentPolygon->second, p1, p2);
        if (p2.contains(p1.back())) {
          using std::swap;
          swap(p1, p2);
        }
        if ((!p1.contains(p2.back()) || id != prevId) && (p1.orientation() > 0) == (hole == 0)) {
          Area& currentSurface = surface[id];
          currentSurface.combine(p1);
          currentSurface.combine(p2);
        } else {
          broken = true;
        }
      }
      if (broken && id == prevId)
        surface.erase(id);
      prevId = id;
    }
    std::unordered_map<int, Polygon> tmp;
    polygon.swap(tmp);
  }

  // Read all zones from the Visum network file.
  std::map<int, Area> zone;
  {
    int id, surfaceId;
    int prevId = INVALID_ID;
    unused(prevId);
    CsvDialect<2> zoneFile(filename + "/BEZIRK.csv");
    zoneFile.read_header(io::ignore_extra_column, "NR", "FLAECHEID");
    while (zoneFile.read_row(id, surfaceId)) {
      assert(id > prevId);
      assert(surface.find(surfaceId) != surface.end());
      zone[id] = surface[surfaceId];
      prevId = id;
    }
  }
  return zone;
}

}
