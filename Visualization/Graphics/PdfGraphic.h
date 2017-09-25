#pragma once

#include <cassert>
#include <string>

#include <cairo/cairo.h>
#include <cairo/cairo-pdf.h>

#include "DataStructures/Geometry/Rectangle.h"
#include "Visualization/Graphic.h"

// A graphic in PDF format. See the description of its base class for further details.
class PdfGraphic : public Graphic {
 public:
  // Constructs a PDF graphic associated with no file.
  PdfGraphic() = default;

  // Constructs a PDF graphic and calls open with the specified arguments.
  PdfGraphic(const std::string& filename, const double width, const double height,
             const Rectangle& userSpace = {{0, 0}, {1, 1}}) {
    open(filename, width, height, userSpace);
  }

  // PDF graphics can only be moved, not copied.
  PdfGraphic(PdfGraphic&& rhs) = default;
  PdfGraphic& operator=(PdfGraphic&& rhs) = default;

  // Opens a PDF graphic file of the specified size (in cm) with the given name.
  void open(const std::string& filename, const double width, const double height,
            const Rectangle& userSpace = {{0, 0}, {1, 1}}) {
    assert(!isOpen());
    deviceWidth = width / 2.54 * 72; // cm to inches to pt
    deviceHeight = height / 2.54 * 72; // cm to inches to pt
    assert(deviceWidth > 0);
    assert(deviceHeight > 0);
    cairo_surface_t* surface;
    surface = cairo_pdf_surface_create((filename + ".pdf").c_str(), deviceWidth, deviceHeight);
    assert(cairo_surface_status(surface) == CAIRO_STATUS_SUCCESS);
    cairoContext = cairo_create(surface);
    assert(cairo_status(cairoContext) == CAIRO_STATUS_SUCCESS);
    cairo_surface_destroy(surface);
    setUserSpace(userSpace);
  }

  // Converts a value measured in points to an equivalent value measured in device-space units.
  virtual double toDeviceSpaceUnits(const double pt) const override {
    return pt;
  }

  // Outputs the current page and inserts a new blank page.
  // Note: This member function does not create a new cairo drawing context.
  virtual void newPage() override {
    assert(isOpen());
    cairo_show_page(cairoContext);
  }
};
