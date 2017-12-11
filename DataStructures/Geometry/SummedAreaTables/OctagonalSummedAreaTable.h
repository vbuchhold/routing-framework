#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "DataStructures/Geometry/Point.h"
#include "DataStructures/Utilities/Matrix.h"

// Implementation of a summed-area table for quickly calculating sums over octagonal areas, as
// described in the expired US Patent 6,507,676.
class OctagonalSummedAreaTable {
 public:
  // Constructs an octagonal summed-area table for the specified matrix.
  explicit OctagonalSummedAreaTable(const Matrix<int>& m)
      : triangular1Sat(m.numRows(), m.numCols()),
        triangular2Sat(m.numRows(), m.numCols()),
        rectangularSat(m.numRows(), m.numCols()) {
    std::vector<int> colSums(m.numCols(), 0);
    for (int y = 0; y != m.numRows(); ++y)
      for (int x = 0; x != m.numCols(); ++x) {
        colSums[x] += m(y, x);
        triangular1Sat(y, x) =
            colSums[x] + (x < m.numCols() - 1 && y > 0 ? triangular1Sat(y - 1, x + 1) : 0);
        triangular2Sat(y, x) = colSums[x] + (x > 0 && y > 0 ? triangular2Sat(y - 1, x - 1) : 0);
        rectangularSat(y, x) = colSums[x] + (x > 0 ? rectangularSat(y, x - 1) : 0);
      }
  }

  // Returns the sum over the octagonal area inscribed in the specified circle.
  int sumOverOctagon(const Point& center, const double radius) const {
    const int x = center.getX();
    const int y = center.getY();
    const int height = std::round(0.92387953251128675613 * radius); // h = r * cos(pi / 8)
    const int side = std::round(0.38268343236508977173 * radius);   // a = r * sin(pi / 8)
    const Point a(x + side, y + height);
    const Point b(x + height + 1, y + side - 1);
    const Point c(x + height, y - side - 1);
    const Point d(x + side - 1, y - height - 2);
    const Point e(x - side + 1, y - height - 2);
    const Point f(x - height, y - side - 1);
    const Point g(x - height - 1, y + side - 1);
    const Point h(x - side, y + height);
    const Point aPrime(x + side - 1, y + height);
    const Point dPrime(x + side - 1, y - height - 1);
    const Point ePrime(x - side, y - height - 1);
    return triangular1SatEntry(a) - triangular1SatEntry(b) -
        triangular2SatEntry(c) + triangular2SatEntry(d) +
        triangular2SatEntry(h) - triangular2SatEntry(g) -
        triangular1SatEntry(f) + triangular1SatEntry(e) +
        rectangularSatEntry(aPrime) - rectangularSatEntry(h) -
        rectangularSatEntry(dPrime) + rectangularSatEntry(ePrime);
  }

 private:
  // Returns the entry at pos of the type 1 triangular SAT.
  int triangular1SatEntry(const Point& pos) const {
    const int x = pos.getX();
    const int y = pos.getY();
    const int numRows = triangular1Sat.numRows();
    const int numCols = triangular1Sat.numCols();
    if (0 <= x && x < numCols && 0 <= y && y < numRows)
      return triangular1Sat(y, x);
    if (x < 0 && -x <= y && y < -x + numRows)
      return triangular1Sat(y + x, 0);
    if (numRows - y <= x && x < numRows - y + numCols - 1 && numRows <= y)
      return triangular1Sat(numRows - 1, x + y - numRows + 1) +
          rectangularSat(numRows - 1, x + y - numRows) -
          (0 < x ? rectangularSat(numRows - 1, x - 1) : 0);
    if (x < numCols && numRows + numCols - x - 1 <= y)
      return rectangularSat(numRows - 1, numCols - 1) -
          (0 < x ? rectangularSat(numRows - 1, x - 1) : 0);
    return 0;
  }

  // Returns the entry at pos of the type 2 triangular SAT.
  int triangular2SatEntry(const Point& pos) const {
    const int x = pos.getX();
    const int y = pos.getY();
    const int numRows = triangular1Sat.numRows();
    const int numCols = triangular1Sat.numCols();
    if (0 <= x && x < numCols && 0 <= y && y < numRows)
      return triangular2Sat(y, x);
    if (numCols <= x && x - numCols < y && y <= x - numCols + numRows)
      return triangular2Sat(y + numCols - x - 1, numCols - 1);
    if (y - numRows < x && x < y - numRows + numCols && numRows <= y)
      return triangular2Sat(numRows - 1, x + numRows - y - 1) -
          rectangularSat(numRows - 1, x + numRows - y - 1) +
          rectangularSat(numRows - 1, std::min(x, numCols - 1));
    if (0 <= x && numRows + x <= y)
      return rectangularSat(numRows - 1, std::min(x, numCols - 1));
    return 0;
  }

  // Returns the entry at pos of the rectangular SAT.
  int rectangularSatEntry(const Point& pos) const {
    const int x = pos.getX();
    const int y = pos.getY();
    const int numRows = triangular1Sat.numRows();
    const int numCols = triangular1Sat.numCols();
    if (0 <= x && 0 <= y)
      return rectangularSat(std::min(y, numRows - 1), std::min(x, numCols - 1));
    return 0;
  }

  Matrix<int> triangular1Sat; // The type 1 triangular summed-area table.
  Matrix<int> triangular2Sat; // The type 2 triangular summed-area table.
  Matrix<int> rectangularSat; // The rectangular summed-area table.
};
