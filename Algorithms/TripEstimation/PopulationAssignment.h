#pragma once

#include <cassert>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <csv.h>

#include "DataStructures/Geometry/CoordinateTransformation.h"
#include "DataStructures/Geometry/Point.h"
#include "DataStructures/Geometry/Rectangle.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Utilities/Matrix.h"
#include "Tools/OpenMP.h"
#include "Tools/Timer.h"

// This class is used to assign a population grid to a graph covering the grid. Each individual in
// a grid cell is assigned to a uniform random vertex in this cell. If the cell does not contain
// vertices, each individual is assigned to a uniform random vertex in the Moore neighborhood of
// range r. The parameter r is increased until the neighborhood is not empty.
template <typename GraphT>
class PopulationAssignment {
 public:
  // Constructs a population assignment procedure with the specified graph and population grid.
  PopulationAssignment(GraphT& graph, const std::string& gridFileName, const int maxRange = 1)
      : graph(graph), gridFileName(gridFileName), maxRange(maxRange) {
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
    const auto primaryCrs = CoordinateTransformation::DHDN_GAUSS_KRUGER_ZONE_3;
    const auto secondaryCrs = CoordinateTransformation::ETRS89_LAEA_EUROPE;
    CoordinateTransformation trans(primaryCrs, secondaryCrs);
    double easting, northing;
    std::vector<Point> cellsByVertex;

    FORALL_VERTICES(graph, v) {
      trans.forward(graph.coordinate(v).getX(), graph.coordinate(v).getY(), easting, northing);
      cellsByVertex.emplace_back(easting / 100, northing / 100);
      boundingBox.extend(cellsByVertex.back());
    }

    // Compute the dimension of the subgrid covered by the graph.
    boundingBox.extend(boundingBox.getSouthWest() - Point(2 * maxRange, 2 * maxRange));
    boundingBox.extend(boundingBox.getNorthEast() + Point(2 * maxRange, 2 * maxRange));
    const auto dim = boundingBox.getNorthEast() - boundingBox.getSouthWest() + Point(1, 1);
    firstVertexInCell.assign(dim.getX() * dim.getY() + 1, 0);
    verticesByCell.assign(graph.numVertices(), 0);
    populationGrid.assign(dim.getY(), dim.getX(), 0);

    FORALL_VERTICES(graph, v) {
      cellsByVertex[v] = cellsByVertex[v] - boundingBox.getSouthWest();
      ++firstVertexInCell[cellId(cellsByVertex[v].getY(), cellsByVertex[v].getX()) + 1];
    }
    std::partial_sum(firstVertexInCell.begin(), firstVertexInCell.end(), firstVertexInCell.begin());

    FORALL_VERTICES(graph, v) {
      const auto id = cellId(cellsByVertex[v].getY(), cellsByVertex[v].getX());
      verticesByCell[firstVertexInCell[id]++] = v;
    }

    for (auto cell = firstVertexInCell.size() - 2; cell > 0; --cell)
      firstVertexInCell[cell] = firstVertexInCell[cell - 1];
    firstVertexInCell[0] = 0;
  }

  // Reads the population grid from file.
  void fillPopulationGrid() {
    Point cell;
    int pop;
    io::CSVReader<3, io::trim_chars<>, io::no_quote_escape<';'>> gridFile(gridFileName);
    gridFile.read_header(io::ignore_extra_column, "x_mp_100m", "y_mp_100m", "Einwohner");
    while (gridFile.read_row(cell.getX(), cell.getY(), pop)) {
      cell.getX() = cell.getX() / 100;
      cell.getY() = cell.getY() / 100;
      if (pop != -1 && boundingBox.contains(cell)) {
        cell = cell - boundingBox.getSouthWest();
        assert(populationGrid(cell.getY(), cell.getX()) == 0);
        populationGrid(cell.getY(), cell.getX()) = pop;
      }
    }
  }

  // Assigns the population in each grid cell to surrounding vertices.
  void assignPopulationToVertices() {
    std::minstd_rand rand(omp_get_thread_num() + 1);
    #pragma omp for collapse(2) nowait
    for (auto i = maxRange; i < populationGrid.numRows() - maxRange; ++i) {
      for (auto j = maxRange; j < populationGrid.numCols() - maxRange; ++j) {
        if (populationGrid(i, j) > 0) {
          auto first = verticesByCell.begin() + firstVertexInCell[cellId(i, j)];
          auto last = verticesByCell.begin() + firstVertexInCell[cellId(i, j) + 1];
          std::vector<int> neighborhood(first, last);

          // Construct the Moore neighborhood of cell (i, j) for increasingly larger ranges.
          for (auto r = 1; neighborhood.empty() && r <= maxRange; ++r) {
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
              #pragma omp atomic
              ++graph.population(neighborhood[dist(rand)]);
            #pragma omp atomic
            assignedPopulation += populationGrid(i, j);
          } else {
            #pragma omp atomic
            unassignedPopulation += populationGrid(i, j);
          }
        }
      }
    }
  }

  // Returns a unique sequential ID for the cell (i, j).
  int cellId(const int i, const int j) const noexcept {
    assert(i >= 0); assert(i < populationGrid.numRows());
    assert(j >= 0); assert(j < populationGrid.numCols());
    return i * populationGrid.numCols() + j;
  }

  GraphT& graph;                  // The graph.
  const std::string gridFileName; // The name of the grid file.
  const int maxRange;             // The upper bound for the range of the Moore neighborhood.

  Rectangle boundingBox;              // A bounding box containing all covered grid cells.
  std::vector<int> firstVertexInCell; // Stores the index of the first vertex in the following list.
  std::vector<int> verticesByCell;    // A list of the vertices ordered by their cell IDs.
  Matrix<int> populationGrid;         // The population grid to be assigned to the graph.

  int assignedPopulation = 0;   // The population that was assigned to a vertex.
  int unassignedPopulation = 0; // The population that could not be assigned to a vertex.
};
