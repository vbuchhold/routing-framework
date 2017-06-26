#pragma once

#include "DataStructures/Geometry/Point.h"
#include "DataStructures/Geometry/Rectangle.h"
#include "Tools/Math.h"

// Returns the normal vector for the line l through p and q, oriented towards the left of l.
inline Point normal(const Point& p, const Point& q) {
  return {p.getY() - q.getY(), q.getX() - p.getX()};
}

// Returns 1, -1, or 0 as r lies to the left of, to the right of, or on the line through p and q.
inline int orientation(const Point& p, const Point& q, const Point& r) {
  return signum(normal(p, q) * (r - p));
}

// Returns true if the line segments pq and rs intersect.
inline bool intersection(const Point& p, const Point& q, const Point& r, const Point& s) {
  const int o1 = orientation(p, q, r);
  const int o2 = orientation(p, q, s);
  if (o1 != o2 && orientation(r, s, p) != orientation(r, s, q))
    return true;
  if (o1 != 0 || o2 != 0)
    return false;
  // All points are collinear.
  Rectangle pq(p);
  Rectangle rs(r);
  pq.extend(q);
  rs.extend(s);
  return pq.intersects(rs);
}
