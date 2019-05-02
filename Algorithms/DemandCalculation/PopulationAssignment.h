#pragma once

#include <cassert>
#include <iostream>
#include <random>
#include <vector>

#include "DataStructures/Containers/BitVector.h"
#include "DataStructures/Geometry/CoordinateTransformation.h"
#include "DataStructures/Geometry/Point.h"
#include "DataStructures/Geometry/Rectangle.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Utilities/Matrix.h"
#include "Tools/LexicalCast.h"
#include "Tools/Math.h"
#include "Tools/OpenMP.h"
#include "Tools/Timer.h"

// This class is used to assign a population grid to a graph covering the grid. Each inhabitant in
// a grid cell is assigned to a uniform random vertex in this cell. If the cell does not contain
// vertices, each inhabitant is assigned to a uniform random vertex in the Moore neighborhood of
// range r. The parameter r is increased until the neighborhood is not empty.
template <typename GraphT, typename GridReaderT>
class PopulationAssignment {
 public:
  // Constructs a population assignment procedure with the specified graph and population grid.
  PopulationAssignment(GraphT& graph, GridReaderT& gridReader, int gridResolution, int maxRange = 1)
      : PopulationAssignment(
            graph, BitVector(graph.numVertices(), true), gridReader, gridResolution, maxRange) {}

  // Constructs a population assignment procedure with the specified graph and population grid.
  PopulationAssignment(GraphT& graph, const BitVector& isVertexInStudyArea,
                       GridReaderT& gridReader, int gridResolution, int maxRange = 1)
      : graph(graph),
        isVertexInStudyArea(isVertexInStudyArea),
        gridReader(gridReader),
        gridResolution(gridResolution),
        maxRange(maxRange) {
    assert(graph.numVertices() == isVertexInStudyArea.size());
    assert(maxRange >= 0);
  }

  void run(const bool verbose = false) {
    Timer timer;
    if (verbose) std::cout << "Assigning vertices to grid cells..." << std::flush;
    assignVerticesToCells();
    if (verbose) std::cout << " done (" << timer.elapsed() << "ms).\n";

    timer.restart();
    if (verbose) std::cout << "Reading population grid from file..." << std::flush;
    fillPopulationGrid();
    if (verbose) std::cout << " done (" << timer.elapsed() << "ms).\n";

    timer.restart();
    if (verbose) std::cout << "Assigning population to vertices..." << std::flush;
    #pragma omp parallel
    assignPopulationToVertices();
    if (verbose) std::cout << " done (" << timer.elapsed() << "ms).\n";

    if (verbose) {
      std::cout << "  Population that was assigned: " << assignedPopulation << "\n";
      std::cout << "  Population that could not be assigned: " << unassignedPopulation << std::endl;
    }
  }

 private:
  // Assigns each vertex to the grid cell containing it.
  void assignVerticesToCells() {
    const auto primaryCrs = CoordinateTransformation::WGS_84;
    const auto secondaryCrs = CoordinateTransformation::ETRS89_LAEA_EUROPE;
    CoordinateTransformation trans(primaryCrs, secondaryCrs);
    double easting, northing;
    std::vector<Point> cellsByVertex;

    for (auto v = isVertexInStudyArea.firstSetBit();
         v != -1;
         v = isVertexInStudyArea.nextSetBit(v)) {
      const auto& latLng = graph.latLng(v);
      trans.forward(toRadians(latLng.lngInDeg()), toRadians(latLng.latInDeg()), easting, northing);
      cellsByVertex.emplace_back(easting / gridResolution, northing / gridResolution);
      boundingBox.extend(cellsByVertex.back());
    }

    // Compute the dimension of the subgrid covered by the graph.
    boundingBox.extend(boundingBox.southWest() - Point(2 * maxRange, 2 * maxRange));
    boundingBox.extend(boundingBox.northEast() + Point(2 * maxRange, 2 * maxRange));
    const auto dim = boundingBox.northEast() - boundingBox.southWest() + Point(1, 1);
    firstVertexInCell.assign(dim.x() * dim.y() + 1, 0);
    verticesByCell.assign(cellsByVertex.size(), 0);
    populationGrid.assign(dim.y(), dim.x(), 0);

    for (auto& cell : cellsByVertex) {
      cell = cell - boundingBox.southWest();
      ++firstVertexInCell[cellId(cell.y(), cell.x()) + 1];
    }
    std::partial_sum(firstVertexInCell.begin(), firstVertexInCell.end(), firstVertexInCell.begin());

    for (auto v = isVertexInStudyArea.firstSetBit(), i = 0;
         v != -1;
         v = isVertexInStudyArea.nextSetBit(v), ++i) {
      const auto id = cellId(cellsByVertex[i].y(), cellsByVertex[i].x());
      verticesByCell[firstVertexInCell[id]++] = v;
    }

    for (auto cell = firstVertexInCell.size() - 2; cell > 0; --cell)
      firstVertexInCell[cell] = firstVertexInCell[cell - 1];
    firstVertexInCell[0] = 0;
  }

  // Reads the population grid from file.
  void fillPopulationGrid() {
    // Set the starting positions of the easting and northing value in an INSPIRE cell code.
    int eastingPos = 0, northingPos = 0;
    switch (gridResolution) {
      case 100:
        eastingPos = 11;
        northingPos = 5;
        break;
      case 1000:
        eastingPos = 9;
        northingPos = 4;
        break;
      default:
        assert(false);
    }

    char* cellCode = nullptr;
    Point cell;
    int pop;

    while (gridReader.read_row(cellCode, pop)) {
      // Decode INSPIRE cell identifier.
      assert(std::strlen(cellCode) > eastingPos);
      assert(cellCode[eastingPos - 1] == 'E');
      assert(cellCode[northingPos - 1] == 'N');
      cellCode[eastingPos - 1] = '\0';
      cell.x() = lexicalCast<int>(cellCode + eastingPos);
      cell.y() = lexicalCast<int>(cellCode + northingPos);

      if (pop != -1 && boundingBox.contains(cell)) {
        cell = cell - boundingBox.southWest();
        populationGrid(cell.y(), cell.x()) += pop;
      }
    }
  }

  // Assigns the population in each grid cell to surrounding vertices.
  void assignPopulationToVertices() {
    std::minstd_rand rand(omp_get_thread_num() + 1);
    std::vector<int> neighborhood;
    std::vector<int> popByVertex(graph.numVertices());
    int assignedPop = 0;
    int unassignedPop = 0;

    #pragma omp for schedule(static, 1024) collapse(2) nowait
    for (auto i = maxRange; i < populationGrid.numRows() - maxRange; ++i) {
      for (auto j = maxRange; j < populationGrid.numCols() - maxRange; ++j) {
        if (populationGrid(i, j) > 0) {
          auto first = verticesByCell.begin() + firstVertexInCell[cellId(i, j)];
          auto last = verticesByCell.begin() + firstVertexInCell[cellId(i, j) + 1];
          neighborhood.assign(first, last);

          // Construct the Moore neighborhood of cell (i, j) for increasingly larger ranges.
          for (auto r = 1; r <= maxRange && neighborhood.empty(); ++r) {
            auto first = verticesByCell.begin() + firstVertexInCell[cellId(i - r, j - r)];
            auto last = verticesByCell.begin() + firstVertexInCell[cellId(i - r, j + r) + 1];
            neighborhood.insert(neighborhood.end(), first, last);
            first = verticesByCell.begin() + firstVertexInCell[cellId(i + r, j - r)];
            last = verticesByCell.begin() + firstVertexInCell[cellId(i + r, j + r) + 1];
            neighborhood.insert(neighborhood.end(), first, last);

            for (auto k = i - r + 1; k <= i + r - 1; ++k) {
              auto first = verticesByCell.begin() + firstVertexInCell[cellId(k, j - r)];
              auto last = verticesByCell.begin() + firstVertexInCell[cellId(k, j - r) + 1];
              neighborhood.insert(neighborhood.end(), first, last);
              first = verticesByCell.begin() + firstVertexInCell[cellId(k, j + r)];
              last = verticesByCell.begin() + firstVertexInCell[cellId(k, j + r) + 1];
              neighborhood.insert(neighborhood.end(), first, last);
            }
          }

          if (!neighborhood.empty()) {
            // Assign each individual in cell (i, j) to a uniform random vertex in the neighborhood.
            std::uniform_int_distribution<> dist(0, neighborhood.size() - 1);
            for (auto k = 0; k < populationGrid(i, j); ++k)
              ++popByVertex[neighborhood[dist(rand)]];
            assignedPop += populationGrid(i, j);
          } else {
            unassignedPop += populationGrid(i, j);
          }
        }
      }
    }

    #pragma omp critical (mergeResults)
    {
      FORALL_VERTICES(graph, v)
        graph.population(v) += popByVertex[v];
      assignedPopulation += assignedPop;
      unassignedPopulation += unassignedPop;
    }
  }

  // Returns a unique sequential ID for the cell (i, j).
  int cellId(const int i, const int j) const noexcept {
    assert(i >= 0); assert(i < populationGrid.numRows());
    assert(j >= 0); assert(j < populationGrid.numCols());
    return i * populationGrid.numCols() + j;
  }

  GraphT& graph;                       // The input graph.
  const BitVector isVertexInStudyArea; // Indicates whether a vertex is in the study area.
  GridReaderT& gridReader;             // A CSV reader for reading the population grid file.
  const int gridResolution;            // The population grid resolution (cell width in meters).
  const int maxRange;                  // The upper bound for the range of the Moore neighborhood.

  Rectangle boundingBox;              // A bounding box containing all covered grid cells.
  std::vector<int> firstVertexInCell; // Stores the index of the first vertex in the following list.
  std::vector<int> verticesByCell;    // A list of the vertices ordered by their cell IDs.
  Matrix<int> populationGrid;         // The population grid to be assigned to the graph.

  int assignedPopulation = 0;   // The population that was assigned to a vertex.
  int unassignedPopulation = 0; // The population that could not be assigned to a vertex.
};
