#pragma once

#include <bitset>
#include <cassert>
#include <limits>

namespace bitwise {

// We add an overloaded version of __builtin_clz.
inline int builtin_clz(const unsigned int x)       { return __builtin_clz(x); }
inline int builtin_clz(const unsigned long x)      { return __builtin_clzl(x); }
inline int builtin_clz(const unsigned long long x) { return __builtin_clzll(x); }

// We add an overloaded version of __builtin_ctz.
inline int builtin_ctz(const unsigned int x)       { return __builtin_ctz(x); }
inline int builtin_ctz(const unsigned long x)      { return __builtin_ctzl(x); }
inline int builtin_ctz(const unsigned long long x) { return __builtin_ctzll(x); }

}

// Returns the bit with the specified index in val.
template <typename UIntT>
inline bool getBit(const UIntT val, const int bitIndex) {
  assert(bitIndex >= 0); assert(bitIndex < std::numeric_limits<UIntT>::digits);
  return (val >> bitIndex) & 1;
}

// Sets the bit with the specified index in val to bitValue.
template <typename UIntT>
inline void setBit(UIntT& val, const int bitIndex, const bool bitValue = true) {
  assert(bitIndex >= 0); assert(bitIndex < std::numeric_limits<UIntT>::digits);
  val ^= (-bitValue ^ val) & (UIntT{1} << bitIndex);
}

// Returns the number of zero-bits preceding the highest-order one-bit in val.
template <typename UIntT>
inline int numLeadingZeros(const UIntT val) {
  // A good compiler will generate a single LZCNT instruction for the following.
  if (val == 0) return std::numeric_limits<UIntT>::digits;
  return bitwise::builtin_clz(val);
}

// Returns the number of zero-bits following the lowest-order one-bit in val.
template <typename UIntT>
inline int numTrailingZeros(const UIntT val) {
  // A good compiler will generate a single TZCNT instruction for the following.
  if (val == 0) return std::numeric_limits<UIntT>::digits;
  return bitwise::builtin_ctz(val);
}

// Returns the number of one-bits in the value representation of val.
template <typename UIntT>
inline int bitCount(const UIntT val) {
  // A good compiler will generate a single POPCNT instruction for the following.
  return std::bitset<std::numeric_limits<UIntT>::digits>(val).count();
}

// Returns the number of one-bits in val that occur before the specified index.
template <typename UIntT>
inline int bitCountBeforeIndex(const UIntT val, const int bitIndex) {
  assert(bitIndex >= 0); assert(bitIndex < std::numeric_limits<UIntT>::digits);
  return bitCount(val & ((UIntT{1} << bitIndex) - 1));
}

// Returns the index of the highest-order one-bit in val, or -1 if val is zero.
template <typename UIntT>
inline int highestOneBit(const UIntT val) {
  return std::numeric_limits<UIntT>::digits - numLeadingZeros(val) - 1;
}

// Returns the index of the lowest-order one-bit in val, or -1 if val is zero.
template <typename UIntT>
inline int lowestOneBit(const UIntT val) {
  if (val == 0) return -1;
  return numTrailingZeros(val);
}

// Returns the index of the highest-order differing bit, or -1 if x and y are equal.
template <typename UIntT>
inline int highestDifferingBit(const UIntT x, const UIntT y) {
  return highestOneBit(x ^ y);
}
