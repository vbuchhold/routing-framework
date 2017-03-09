#pragma once

#include <cstdint>

// Statistics about an iterative all-or-nothing assignment, including checksums and running times.
struct AllOrNothingAssignmentStats {
  // Constructs a struct collecting statistics about an iterative all-or-nothing assignment.
  AllOrNothingAssignmentStats()
      : lastChecksum(0),
        totalChecksum(0),
        lastCustomizationTime(0),
        lastQueryTime(0),
        totalPreprocessingTime(0),
        totalCustomizationTime(0),
        totalQueryTime(0),
        numIterations(0) {}

  // Adds the values from the last iteration to the totals.
  void addLastValuesToTotals() {
    totalChecksum += lastChecksum;
    totalCustomizationTime += lastCustomizationTime;
    totalQueryTime += lastQueryTime;
  }

  int64_t lastChecksum;  // The sum of the distances computed in the last iteration.
  int64_t totalChecksum; // The total sum of distances computed.

  int lastCustomizationTime; // The time spent on customization in the last iteration.
  int lastQueryTime;         // The time spent on queries in the last iteration.

  int totalPreprocessingTime; // The total time spent on preprocessing.
  int totalCustomizationTime; // The total time spent on customization.
  int totalQueryTime;         // The total time spent on queries.

  int numIterations; // The number of iterations performed.
};
