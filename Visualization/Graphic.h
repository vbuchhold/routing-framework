#pragma once

#include <cassert>

#include <cairo/cairo.h>

#include "DataStructures/Geometry/Rectangle.h"

// A graphic on which geometric objects can be drawn. To position objects on the graphic independent
// of its type and size, each graphic maintains a user-space coordinate system. The Cartesian user
// space is specified by its smallest and largest x- and y-coordinate. The point with the smallest
// (largest) coordinates is transformed to the lower left (upper right) corner of the graphic.
// However, the user space is not distorted, i.e., one user-space unit in x-direction is in device
// space always as long as a user-space unit in y-direction. This is achieved by adding horizontal
// or vertical padding if needed. This is an abstract base class. Derived classes provide graphics
// of concrete types.
class Graphic {
 public:
  // Constructs a graphic associated with no file.
  Graphic() : cairoContext(nullptr), deviceWidth(0), deviceHeight(0) {}

  // Outputs the current page (if any) and releases resources.
  virtual ~Graphic() {
    cairo_destroy(cairoContext);
  }

  // Returns true if this graphic is associated with a file and ready to be drawn on.
  bool isOpen() const {
    return cairoContext != nullptr;
  }

  // Outputs the current page (if any) and releases resources.
  virtual void close() {
    cairo_destroy(cairoContext);
    cairoContext = nullptr;
  }

  // Sets the user space according to the specified rectangle.
  void setUserSpace(const Rectangle& uspace) {
    assert(isOpen());
    const double uspaceWidth = uspace.getNorthEast().getX() - uspace.getSouthWest().getX();
    const double uspaceHeight = uspace.getNorthEast().getY() - uspace.getSouthWest().getY();
    assert(uspaceWidth > 0);
    assert(uspaceHeight > 0);
    double scale, offsetX = 0, offsetY = 0;
    if (deviceWidth * uspaceHeight <= deviceHeight * uspaceWidth) {
      scale = deviceWidth / uspaceWidth;
      offsetY = (deviceHeight - scale * uspaceHeight) / 2.0;
    } else {
      scale = deviceHeight / uspaceHeight;
      offsetX = (deviceWidth - scale * uspaceWidth) / 2.0;
    }
    cairo_identity_matrix(cairoContext);
    cairo_translate(cairoContext, offsetX, offsetY);
    cairo_scale(cairoContext, scale, scale);
    cairo_translate(cairoContext, -uspace.getSouthWest().getX(), uspace.getNorthEast().getY());
    cairo_scale(cairoContext, 1, -1);
    userSpace = uspace;
  }

  // Returns the user space.
  Rectangle getUserSpace() const {
    return userSpace;
  }

  // Returns the encapsulated cairo drawing context.
  cairo_t* getCairoContext() {
    return cairoContext;
  }

  // Converts a value measured in points to an equivalent value measured in device-space units.
  virtual double toDeviceSpaceUnits(const double pt) const = 0;

  // Converts a value measured in points to an equivalent value measured in user-space units.
  double toUserSpaceUnits(const double pt) {
    assert(isOpen());
    double x = toDeviceSpaceUnits(pt), y = 0;
    cairo_device_to_user_distance(cairoContext, &x, &y);
    return x;
  }

  // Outputs the current page and inserts a new blank page.
  // Note: This member function may or may not create a new cairo drawing context.
  virtual void newPage() = 0;

 protected:
  cairo_t* cairoContext; // A cairo drawing context used to draw on the graphic.
  Rectangle userSpace;   // The user space specified as a rectangle.
  double deviceWidth;    // The width of the graphic in points or pixels, depending on the type.
  double deviceHeight;   // The height of the graphic in points or pixels, depending on the type.
};
