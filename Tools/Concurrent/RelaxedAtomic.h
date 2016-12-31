#pragma once

#include <atomic>

// A wrapper around std::atomic. It has the same interface as std::atomic. The only difference is
// that this class performs atomic operations with relaxed ordering by default, whereas std::atomic
// performs them with sequentially consistent ordering by default.
template <typename T>
class RelaxedAtomic : public std::atomic<T> {
 public:
  void store(T desired, std::memory_order order = std::memory_order_relaxed) noexcept {
    std::atomic<T>::store(desired, order);
  }

  T load(std::memory_order order = std::memory_order_relaxed) const noexcept {
    return std::atomic<T>::load(order);
  }

  operator T() const noexcept {
    return load();
  }

  T operator=(T desired) noexcept {
    store(desired);
    return desired;
  }
};
