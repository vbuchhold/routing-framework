#pragma once

#include <chrono>

// A timer to measure how long some code takes to execute.
class Timer {
 public:
  // Constructs a timer and starts it.
  Timer() : startTime(std::chrono::steady_clock::now()) {}

  // Returns the time elapsed since the timer was started.
  template <typename UnitT = std::chrono::milliseconds>
  int elapsed() const {
    const auto now = std::chrono::steady_clock::now();
    return std::chrono::duration_cast<UnitT>(now - startTime).count();
  }

  // Restarts the timer.
  void restart() {
    startTime = std::chrono::steady_clock::now();
  }

 private:
  std::chrono::steady_clock::time_point startTime; // Time point when the timer was started.
};
