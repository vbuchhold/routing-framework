#pragma once

#include <vectorclass/vectorclass.h>

// Returns the horizontal minimum amongst the packed 32-bit integers in the specified vector.
int horizontal_min(Vec4i packedMin) {
  packedMin = min(packedMin, permute4i<2, 3, -256, -256>(packedMin));
  packedMin = min(packedMin, permute4i<1, -256, -256, -256>(packedMin));
  return packedMin.extract(0);
}

// Returns the horizontal maximum amongst the packed 32-bit integers in the specified vector.
int horizontal_max(Vec4i packedMax) {
  packedMax = max(packedMax, permute4i<2, 3, -256, -256>(packedMax));
  packedMax = max(packedMax, permute4i<1, -256, -256, -256>(packedMax));
  return packedMax.extract(0);
}

// Returns the horizontal minimum amongst the packed 32-bit integers in the specified vector.
int horizontal_min(Vec8i packedMin) {
  packedMin = min(packedMin, permute8i<4, 5, 6, 7, -256, -256, -256, -256>(packedMin));
  packedMin = min(packedMin, permute8i<2, 3, -256, -256, -256, -256, -256, -256>(packedMin));
  packedMin = min(packedMin, permute8i<1, -256, -256, -256, -256, -256, -256, -256>(packedMin));
  return packedMin.extract(0);
}

// Returns the horizontal maximum amongst the packed 32-bit integers in the specified vector.
int horizontal_max(Vec8i packedMax) {
  packedMax = max(packedMax, permute8i<4, 5, 6, 7, -256, -256, -256, -256>(packedMax));
  packedMax = max(packedMax, permute8i<2, 3, -256, -256, -256, -256, -256, -256>(packedMax));
  packedMax = max(packedMax, permute8i<1, -256, -256, -256, -256, -256, -256, -256>(packedMax));
  return packedMax.extract(0);
}
