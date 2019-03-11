#pragma once

#include <limits>

// Statistics about a Frank-Wolfe assignment, including times and measures of solution quality.
struct FrankWolfeAssignmentStats {
  // Constructs a struct collecting statistics about a Frank-Wolfe assignment.
  FrankWolfeAssignmentStats()
      : prevTotalTraversalCost(0),
        prevTotalPathCost(0),
        prevRelGap(std::numeric_limits<double>::infinity()),
        lastLineSearchTime(0),
        lastRunningTime(0),
        totalLineSearchTime(0),
        totalRunningTime(0) {}

  // Resets the values from the last iteration.
  void startIteration() {
    prevTotalTraversalCost = 0;
    prevTotalPathCost = 0;
  }

  // Adds the values from the last iteration to the totals.
  void finishIteration() {
    totalLineSearchTime += lastLineSearchTime;
    totalRunningTime += lastRunningTime;
  }

  double prevTotalTraversalCost; // The total traversal cost after the previous iteration.
  double prevTotalPathCost;      // The total path cost after the previous iteration.
  double prevRelGap;             // The relative gap after the previous iteration.

  int lastLineSearchTime; // The time spent on the line search in the last iteration.
  int lastRunningTime;    // The running time for the last iteration.

  int totalLineSearchTime; // The total time spent on the line search.
  int totalRunningTime;    // The total running time.
};
