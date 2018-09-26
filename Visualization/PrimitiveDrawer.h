#pragma once

#include <cassert>
#include <string>

#include <cairo/cairo.h>

#include "DataStructures/Geometry/Point.h"
#include "DataStructures/Geometry/Polygon.h"
#include "Visualization/Color.h"
#include "Visualization/Graphic.h"

// Some predefined line width values.
struct LineWidth {
  static constexpr double ULTRA_THIN = 0.1;
  static constexpr double VERY_THIN = 0.2;
  static constexpr double THIN = 0.4;
  static constexpr double SEMITHICK = 0.6;
  static constexpr double THICK = 0.8;
  static constexpr double VERY_THICK = 1.2;
  static constexpr double ULTRA_THICK = 1.6;
};

// A facility for drawing geometric primitives like lines and polygons on a graphic.
class PrimitiveDrawer {
 public:
  // Constructs a drawer for geometric primitives.
  PrimitiveDrawer(Graphic* const graphic = nullptr,
                  const Color& color = KIT_BLACK, const double width = LineWidth::THIN)
      : currentGraphic(graphic), currentColor(color), currentLineWidth(width) {
    initCairoContext();
  }

  // Initializes the cairo drawing context of the current graphic.
  void initCairoContext() {
    if (currentGraphic != nullptr && currentGraphic->isOpen()) {
      cairo_t* const ctx = currentGraphic->getCairoContext();
      setColor(currentColor);
      setLineWidth(currentLineWidth);
      cairo_set_line_cap(ctx, CAIRO_LINE_CAP_ROUND);
      cairo_set_font_size(ctx, currentGraphic->toUserSpaceUnits(10));
    }
  }

  // Sets the graphic to be drawn on.
  void setGraphic(Graphic& graphic) {
    currentGraphic = &graphic;
    initCairoContext();
  }

  // Sets the color to be used for drawing.
  void setColor(const Color& color) {
    currentColor = color;
    if (currentGraphic != nullptr && currentGraphic->getCairoContext() != nullptr)
      cairo_set_source_rgba(
          currentGraphic->getCairoContext(),
          currentColor.normalizedRed(), currentColor.normalizedGreen(),
          currentColor.normalizedBlue(), currentColor.normalizedAlpha());
  }

  // Sets the line width (in points) to be used for drawing.
  void setLineWidth(const double width) {
    currentLineWidth = width;
    if (currentGraphic != nullptr && currentGraphic->getCairoContext() != nullptr) {
      const double userSpaceUnits = currentGraphic->toUserSpaceUnits(currentLineWidth);
      cairo_set_line_width(currentGraphic->getCairoContext(), userSpaceUnits);
    }
  }

  // Draws a line between the specified points, possibly setting the current color and line width.
  void drawLine(const Point& from, const Point& to) {
    assert(currentGraphic != nullptr); assert(currentGraphic->isOpen());
    cairo_t* const ctx = currentGraphic->getCairoContext();
    cairo_move_to(ctx, from.getX(), from.getY());
    cairo_line_to(ctx, to.getX(), to.getY());
    cairo_stroke(ctx);
  }
  void drawLine(const Point& from, const Point& to, const Color& color) {
    setColor(color);
    drawLine(from, to);
  }
  void drawLine(const Point& from, const Point& to, const double width) {
    setLineWidth(width);
    drawLine(from, to);
  }
  void drawLine(const Point& from, const Point& to, const Color& color, const double width) {
    setColor(color);
    setLineWidth(width);
    drawLine(from, to);
  }

  // Draws the specified polygon, possibly setting the current color and line width.
  void drawPolygon(const Polygon& polygon) {
    assert(currentGraphic != nullptr); assert(currentGraphic->isOpen());
    if (polygon.empty())
      return;
    cairo_t* const ctx = currentGraphic->getCairoContext();
    cairo_move_to(ctx, polygon[0].getX(), polygon[0].getY());
    for (int i = 1; i < polygon.size(); ++i)
      cairo_line_to(ctx, polygon[i].getX(), polygon[i].getY());
    cairo_close_path(ctx);
    cairo_stroke(ctx);
  }
  void drawPolygon(const Polygon& polygon, const Color& color) {
    setColor(color);
    drawPolygon(polygon);
  }
  void drawPolygon(const Polygon& polygon, const double width) {
    setLineWidth(width);
    drawPolygon(polygon);
  }
  void drawPolygon(const Polygon& polygon, const Color& color, const double width) {
    setColor(color);
    setLineWidth(width);
    drawPolygon(polygon);
  }

  // Fills a circle with the specified center and radius, possibly setting the current color.
  void fillCircle(const Point& center, const double radius) {
    assert(currentGraphic != nullptr); assert(currentGraphic->isOpen());
    cairo_t* const ctx = currentGraphic->getCairoContext();
    cairo_arc(ctx, center.getX(), center.getY(), radius, 0, 2 * PI);
    cairo_fill(ctx);
  }
  void fillCircle(const Point& center, const double radius, const Color& color) {
    setColor(color);
    fillCircle(center, radius);
  }

  // Draws the text given by the specified string centered on the specified position.
  void drawString(const std::string& str, const Point& pos) {
    assert(currentGraphic != nullptr); assert(currentGraphic->isOpen());
    cairo_t* const ctx = currentGraphic->getCairoContext();
    cairo_text_extents_t strMetrics;
    cairo_text_extents(ctx, str.c_str(), &strMetrics);
    const auto x = pos.getX() - strMetrics.width / 2 - strMetrics.x_bearing;
    const auto y = pos.getY() + strMetrics.height / 2 + strMetrics.y_bearing;
    cairo_move_to(ctx, x, y);
    cairo_scale(ctx, 1, -1);
    cairo_show_text(ctx, str.c_str());
    cairo_scale(ctx, 1, -1);
  }

 private:
  Graphic* currentGraphic; // The graphic to be drawn on.
  Color currentColor;      // The color currently used for drawing.
  double currentLineWidth; // The line width currently used for drawing.
};
