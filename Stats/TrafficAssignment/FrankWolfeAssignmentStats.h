#pragma once

// Statistics about a Frank-Wolfe assignment, including times and measures of solution quality.
struct FrankWolfeAssignmentStats {
  // Constructs a struct collecting statistics about a Frank-Wolfe assignment.
  FrankWolfeAssignmentStats()
      : changeInEdgeFlows(0),
        totalTravelCost(0),
        lastRunningTime(0),
        totalRunningTime(0) {}

  // Resets the values from the last iteration.
  void startIteration() {
    changeInEdgeFlows = 0;
    totalTravelCost = 0;
  }

  // Adds the values from the last iteration to the totals.
  void finishIteration() {
    totalRunningTime += lastRunningTime;
  }

  float changeInEdgeFlows; // The change in the edge flows between the last two iterations.
  float totalTravelCost;   // The total travel cost resulting from the current edge flows.

  int lastRunningTime;  // The running time for the last iteration.
  int totalRunningTime; // The total running time.
};
