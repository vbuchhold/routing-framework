#pragma once

#include <cassert>
#include <iostream>
#include <ostream>

#include <omp.h>

// A textual indicator of progress towards some goal.
class ProgressBar {
 public:
  // Constructs an uninitialized progress bar.
  explicit ProgressBar(const bool verbose = true, std::ostream& os = std::cout)
      : os(os), numSteps(0), stepsDone(0), percentageDone(0),
        percentageOutputInterval(20), dotOutputInterval(5), verbose(verbose) {}

  // Constructs a progress bar with the specified number of steps.
  explicit ProgressBar(const int numSteps, const bool verbose = true, std::ostream& os = std::cout)
      : ProgressBar(verbose, os) {
    init(numSteps);
  }

  // Initialize the progress bar with the specified number of steps.
  void init(const int steps) {
    assert(steps >= 0);
    numSteps = steps;
    stepsDone = 0;
    percentageDone = 0;
    if (verbose)
      os << "0% " << std::flush;
  }

  // Set the percentage points between two printed percentages.
  void setPercentageOutputInterval(const int points) {
    percentageOutputInterval = points;
  }

  // Set the percentage points between two printed dots.
  void setDotOutputInterval(const int points) {
    dotOutputInterval = points;
  }

  // Advances the progress bar to the specified step.
  void advanceTo(const int step) {
    if (!verbose)
      return;
    assert(step >= stepsDone); assert(step <= numSteps);
    stepsDone = step;
    print(stepsDone * 100l / numSteps);
  }

  // Advances the progress bar to 100 %.
  void finish() {
    advanceTo(numSteps);
  }

  // Advances the progress bar by one step.
  void operator++() {
    if (!verbose)
      return;
    int done;
    #pragma omp atomic capture
    done = ++stepsDone;
#ifdef _OPENMP
    if (omp_get_thread_num() == 0)
#endif
      print(done * 100l / numSteps);
  }

  // Advances the progress bar by the specified number of steps.
  void operator+=(const int steps) {
    assert(steps >= 0);
    if (!verbose)
      return;
    int done;
    #pragma omp atomic capture
    done = stepsDone += steps;
#ifdef _OPENMP
    if (omp_get_thread_num() == 0)
#endif
      print(done * 100l / numSteps);
  }

 private:
  // Prints the progress bar until the specified percentage.
  void print(const int until) {
    assert(until <= 100);
    for (int i = percentageDone + 1; i <= until; ++i)
      if (i % percentageOutputInterval == 0)
        os << " " << i << "% " << std::flush;
      else if (i % dotOutputInterval == 0)
        os << "." << std::flush;
    percentageDone = until;
  }

  std::ostream& os; // The output stream the progress bar is printed to.

  int numSteps;       // The number of steps that have to be done.
  int stepsDone;      // The number of steps that have already been done.
  int percentageDone; // The percentage that has already been done.

  int percentageOutputInterval; // Percentage points between two printed percentages.
  int dotOutputInterval;        // Percentage points between two printed dots.
  bool verbose;                 // Indicates if the progress bar should be printed.
};
