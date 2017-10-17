#pragma once

#include <cassert>
#include <limits>

namespace bitwise {

#if defined __GNUC__
// Returns the number of leading 0-bits in x. If x is zero, the result is undefined.
inline int clz(const unsigned int x) { return __builtin_clz(x); }
inline int clz(const unsigned long x) { return __builtin_clzl(x); }
inline int clz(const unsigned long long x) { return __builtin_clzll(x); }
#endif

}

// Returns the i-th bit in x.
template <typename T>
inline bool getBit(const T x, const int i) {
  assert(i >= 0); assert(i < std::numeric_limits<T>::digits);
  return (x >> i) & 1;
}

// Sets the i-th bit in x to val.
template <typename T>
inline void setBit(T& x, const int i, const bool val = true) {
  assert(i >= 0); assert(i < std::numeric_limits<T>::digits);
  x ^= (-val ^ x) & (static_cast<T>(1) << i);
}

// Returns the index of the most significant 1-bit in x, or -1 if x is zero.
template <typename T>
inline int mostSignificantOneBit(T x) {
#if defined __GNUC__
  return x ? std::numeric_limits<T>::digits - bitwise::clz(x) - 1 : -1;
#else
  int idx = -1;
  while (x) {
    x >>= 1;
    ++idx;
  }
  return idx;
#endif
}

// Returns the index of the most significant differing bit, or -1 if x and y are the same.
template <typename T>
inline int mostSignificantDifferingBit(const T x, const T y) {
  return mostSignificantOneBit(x ^ y);
}
