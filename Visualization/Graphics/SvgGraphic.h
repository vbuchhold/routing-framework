#pragma once

#include <cassert>
#include <string>

#include <cairo/cairo.h>
#include <cairo/cairo-svg.h>

#include "DataStructures/Geometry/Rectangle.h"
#include "Visualization/Graphic.h"

// A graphic in SVG format. See the description of its base class for further details.
class SvgGraphic : public Graphic {
 public:
  // Constructs a SVG graphic associated with no file.
  SvgGraphic() : page(0) {}

  // Constructs a SVG graphic and calls open with the specified arguments.
  SvgGraphic(const std::string& filename, const double width, const double height,
             const Rectangle& userSpace = {{0, 0}, {1, 1}}) {
    open(filename, width, height, userSpace);
  }

  // SVG graphics can only be moved, not copied.
  SvgGraphic(SvgGraphic&& rhs) = default;
  SvgGraphic& operator=(SvgGraphic&& rhs) = default;

  // Opens a SVG graphic file of the specified size (in cm) with the given name.
  void open(const std::string& filename, const double width, const double height,
            const Rectangle& userSpace = {{0, 0}, {1, 1}}) {
    assert(!isOpen());
    deviceWidth = width / 2.54 * 72; // cm to inches to pt
    deviceHeight = height / 2.54 * 72; // cm to inches to pt
    assert(deviceWidth > 0);
    assert(deviceHeight > 0);
    cairo_surface_t* surface;
    surface = cairo_svg_surface_create((filename + ".svg").c_str(), deviceWidth, deviceHeight);
    assert(cairo_surface_status(surface) == CAIRO_STATUS_SUCCESS);
    cairoContext = cairo_create(surface);
    assert(cairo_status(cairoContext) == CAIRO_STATUS_SUCCESS);
    cairo_surface_destroy(surface);
    setUserSpace(userSpace);
    baseFilename = filename;
    page = 1;
  }

  // Converts a value measured in points to an equivalent value measured in device-space units.
  virtual double toDeviceSpaceUnits(const double pt) const override {
    return pt;
  }

  // Outputs the current page and inserts a new blank page.
  // Note: This member function creates a new cairo drawing context.
  virtual void newPage() override {
    assert(isOpen());
    close();
    const std::string name = baseFilename + "_" + std::to_string(++page) + ".svg";
    cairo_surface_t* surface = cairo_svg_surface_create(name.c_str(), deviceWidth, deviceHeight);
    assert(cairo_surface_status(surface) == CAIRO_STATUS_SUCCESS);
    cairoContext = cairo_create(surface);
    assert(cairo_status(cairoContext) == CAIRO_STATUS_SUCCESS);
    cairo_surface_destroy(surface);
    setUserSpace(userSpace);
  }

 private:
  std::string baseFilename; // The base name of the associated graphic file.
  int page;                 // The number of the current page.
};
