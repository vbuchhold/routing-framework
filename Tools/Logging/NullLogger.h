#pragma once

// A logger that discards the data written to it, avoiding any overhead at runtime.
class NullLogger {
 public:
  template <typename T>
  explicit NullLogger(const T&) noexcept {}

  explicit operator bool() const noexcept {
    return true;
  }

  template <typename T>
  NullLogger& operator<<(const T&) noexcept {
    return *this;
  }
};
