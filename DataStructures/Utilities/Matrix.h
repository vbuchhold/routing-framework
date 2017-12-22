#pragma once

#include <cassert>
#include <vector>

// A rectangular array of elements in rows and columns.
template <typename E>
class Matrix {
 public:
  // Constructs a matrix with the specified number of rows and columns.
  Matrix(const int numRows, const int numCols, const E& val = E())
      : data(numRows * numCols, val), rowCount(numRows), colCount(numCols) {
    assert(numRows >= 0);
    assert(numCols >= 0);
  }

  // Returns the element of the matrix in row i and column j.
  E& operator()(const int i, const int j) {
    assert(i >= 0); assert(i < rowCount);
    assert(j >= 0); assert(j < colCount);
    return data[i * colCount + j];
  }

  // Returns the element of the matrix in row i and column j.
  const E& operator()(const int i, const int j) const {
    assert(i >= 0); assert(i < rowCount);
    assert(j >= 0); assert(j < colCount);
    return data[i * colCount + j];
  }

  // Returns the number of rows in the matrix.
  int numRows() const {
    return rowCount;
  }

  // Returns the number of columns in the matrix.
  int numCols() const {
    return colCount;
  }

 private:
  std::vector<E> data; // The elements of the matrix in row-major order.
  const int rowCount;  // The number of rows in the matrix.
  const int colCount;  // The number of columns in the matrix.
};
