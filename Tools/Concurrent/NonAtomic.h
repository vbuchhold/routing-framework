#pragma once

#include <atomic>

// This class represents a plain old non-atomic variable of type T. It has the same interface as
// std::atomic and allows non-atomics to be used where atomics are expected.
template <typename T>
class NonAtomic {
 public:
  void store(T desired, std::memory_order /*order*/ = std::memory_order_seq_cst) noexcept {
    value = desired;
  }

  T load(std::memory_order /*order*/ = std::memory_order_seq_cst) const noexcept {
    return value;
  }

  operator T() const noexcept {
    return load();
  }

  T operator=(T desired) noexcept {
    store(desired);
    return desired;
  }

 private:
  T value; // The value of the variable.
};
