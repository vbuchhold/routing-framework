#pragma once

#include <functional>

#include "DataStructures/Geometry/Point.h"
#include "DataStructures/Graph/Graph.h"
#include "Visualization/Color.h"
#include "Visualization/Graphic.h"
#include "Visualization/PrimitiveDrawer.h"

// A facility for drawing graphs. Each edge is drawn as a straight line segment. The color (width)
// of an edge is determined by a user-provided callable object that takes an edge ID as argument
// and returns the corresponding color (width).
class GraphDrawer : public PrimitiveDrawer {
 public:
  // Constructs a drawer for graphs.
  explicit GraphDrawer(
      Graphic* const graphic = nullptr,
      const std::function<Point(int)>& coordinate = {},
      const std::function<Color(int)>& color = [](const int /*e*/) { return KIT_BLACK; },
      const std::function<double(int)>& width = [](const int /*e*/) { return LineWidth::THIN; })
      : PrimitiveDrawer(graphic), getCoordinate(coordinate), getColor(color), getLineWidth(width) {}

  // Sets the callable object determining the coordinate of a vertex.
  void setCoordinate(const std::function<Point(int)>& coordinate) {
    getCoordinate = coordinate;
  }

  // Sets the callable object determining the color of an edge.
  void setColor(const std::function<Color(int)>& color) {
    getColor = color;
  }

  // Sets the callable object determining the width of an edge.
  void setLineWidth(const std::function<double(int)>& width) {
    getLineWidth = width;
  }

  // Draws the specified graph.
  template <typename GraphT>
  void drawGraph(const GraphT& graph) {
    FORALL_VALID_EDGES(graph, u, e)
      drawLine(getCoordinate(u), getCoordinate(graph.edgeHead(e)), getColor(e), getLineWidth(e));
  }

 private:
  std::function<Point(int)> getCoordinate; // The function determining the coordinate of a vertex.
  std::function<Color(int)> getColor;      // The function determining the color of an edge.
  std::function<double(int)> getLineWidth; // The function determining the width of an edge.
};
