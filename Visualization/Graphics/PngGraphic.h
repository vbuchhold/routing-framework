#pragma once

#include <cassert>
#include <cmath>
#include <string>

#include <cairo/cairo.h>

#include "DataStructures/Geometry/Rectangle.h"
#include "Visualization/Graphic.h"

// A graphic in PNG format with 381dpi. See the description of its base class for further details.
class PngGraphic : public Graphic {
 public:
  // Constructs a PNG graphic associated with no file.
  PngGraphic() : page(0) {}

  // Constructs a PNG graphic and calls open with the specified arguments.
  PngGraphic(const std::string& filename, const double width, const double height,
             const Rectangle& userSpace = {{0, 0}, {1, 1}}) {
    open(filename, width, height, userSpace);
  }

  // PNG graphics can only be moved, not copied.
  PngGraphic(PngGraphic&& rhs) = default;
  PngGraphic& operator=(PngGraphic&& rhs) = default;

  // Outputs the current page (if any) and releases resources.
  virtual ~PngGraphic() override {
    if (isOpen()) {
      const std::string name = baseFilename + (page > 1 ? "_" + std::to_string(page) : "");
      cairo_surface_write_to_png(cairo_get_target(cairoContext), (name + ".png").c_str());
    }
  }

  // Opens a PNG graphic file of the specified size (in cm) with the given name.
  void open(const std::string& filename, const double width, const double height,
            const Rectangle& userSpace = {{0, 0}, {1, 1}}) {
    assert(!isOpen());
    deviceWidth = std::round(width / 2.54 * 381); // cm to inches to pixels
    deviceHeight = std::round(height / 2.54 * 381); // cm to inches to pixels
    assert(deviceWidth > 0);
    assert(deviceHeight > 0);
    cairo_surface_t* surface;
    surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, deviceWidth, deviceHeight);
    assert(cairo_surface_status(surface) == CAIRO_STATUS_SUCCESS);
    cairoContext = cairo_create(surface);
    assert(cairo_status(cairoContext) == CAIRO_STATUS_SUCCESS);
    cairo_surface_destroy(surface);
    setUserSpace(userSpace);
    baseFilename = filename;
    page = 1;
  }

  // Outputs the current page (if any) and releases resources.
  virtual void close() override {
    if (isOpen()) {
      const std::string name = baseFilename + (page > 1 ? "_" + std::to_string(page) : "");
      cairo_surface_write_to_png(cairo_get_target(cairoContext), (name + ".png").c_str());
      Graphic::close();
    }
  }

  // Converts a value measured in points to an equivalent value measured in device-space units.
  virtual double toDeviceSpaceUnits(const double pt) const override {
    return pt / 72 * 381; // pt to inches to pixels
  }

  // Outputs the current page and inserts a new blank page.
  // Note: This member function creates a new cairo drawing context.
  virtual void newPage() override {
    assert(isOpen());
    close();
    cairo_surface_t* surface;
    surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, deviceWidth, deviceHeight);
    assert(cairo_surface_status(surface) == CAIRO_STATUS_SUCCESS);
    cairoContext = cairo_create(surface);
    assert(cairo_status(cairoContext) == CAIRO_STATUS_SUCCESS);
    cairo_surface_destroy(surface);
    setUserSpace(userSpace);
    ++page;
  }

 private:
  std::string baseFilename; // The base name of the associated graphic file.
  int page;                 // The number of the current page.
};
