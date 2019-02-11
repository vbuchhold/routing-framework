#pragma once

#include <cassert>
#include <string>

#include <proj.h>

// This class provides coordinate transformation between two coordinate reference systems, which
// are specified as EPSG codes.
class CoordinateTransformation {
 public:
  // EPSG codes for some widely used coordinate systems.
  static constexpr int ETRS89_LAEA_EUROPE = 3035;
  static constexpr int WGS_84 = 4326;
  static constexpr int OSGB_1936_BRITISH_NATIONAL_GRID = 27700;
  static constexpr int DHDN_GAUSS_KRUGER_ZONE_3 = 31467;

  // Constructs a transformation between the two coordinate systems specified as EPSG codes.
  CoordinateTransformation(const int primaryCrs, const int secondaryCrs) {
    const auto projStr = "+proj=pipeline +step +inv +init=epsg:" + std::to_string(primaryCrs) +
        " +step +init=epsg:" + std::to_string(secondaryCrs);
    threadingContext = proj_context_create();
    transformation = proj_create(threadingContext, projStr.c_str());
    assert(transformation != 0);
  }

  // Copy constructor.
  CoordinateTransformation(const CoordinateTransformation& other) {
    threadingContext = proj_context_create();
    transformation = proj_create(threadingContext, proj_pj_info(other.transformation).definition);
    assert(transformation != 0);
  }

  // Releases all resources.
  ~CoordinateTransformation() {
    proj_destroy(transformation);
    proj_context_destroy(threadingContext);
  }

  // Transforms a coordinate given in the primary CRS into a coordinate in the secondary CRS.
  void forward(const double srcX, const double srcY, double& dstX, double& dstY) {
    const auto result = proj_trans(transformation, PJ_FWD, {{srcX, srcY, 0, 0}});
    dstX = result.xy.x;
    dstY = result.xy.y;
  }

  // Transforms a coordinate given in the secondary CRS into a coordinate in the primary CRS.
  void reverse(const double srcX, const double srcY, double& dstX, double& dstY) {
    const auto result = proj_trans(transformation, PJ_INV, {{srcX, srcY, 0, 0}});
    dstX = result.xy.x;
    dstY = result.xy.y;
  }

 private:
  PJ_CONTEXT* threadingContext; // The PROJ threading context.
  PJ* transformation;           // The PROJ transformation object.
};
