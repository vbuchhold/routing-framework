#pragma once

// Atomically replaces the value of the first argument with the smaller of two values.
template <typename T>
inline void atomicFetchMin(T& a, const T& b) {
  T expect = a;
  while (b < expect &&
         !__atomic_compare_exchange_n(&a, &expect, b, true, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) {}
}
