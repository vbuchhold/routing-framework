#pragma once

#include <cassert>
#include <algorithm>
#include <cstring>
#include <string>

namespace str {

// Returns an iterator referring to the first character in the specified string.
inline std::string::iterator begin(std::string& string) { return string.begin(); }
inline std::string::const_iterator begin(const std::string& string) { return string.begin(); }
inline char* begin(char* string) { return string; }
inline const char* begin(const char* string) { return string; }

// Returns the length of the specified string.
inline int length(const std::string& string) { return string.length(); }
inline int length(const char* const string) { return std::strlen(string); }

// Returns a pointer referring to the first character in the specified string.
inline const char* cStr(const std::string& string) { return string.c_str(); }
inline const char* cStr(const char* string) { return string; }

}

// Tests if the specified string starts with the specified prefix.
template <typename StringT1, typename StringT2>
inline bool startsWith(const StringT1& string, const StringT2& prefix) {
  const int strlen = str::length(string);
  const int prelen = str::length(prefix);
  if (strlen < prelen)
    return false;
  return std::equal(str::begin(string), str::begin(string) + prelen, str::begin(prefix));
}

// Tests if the specified string ends with the specified suffix.
template <typename StringT1, typename StringT2>
inline bool endsWith(const StringT1& string, const StringT2& suffix) {
  const int strlen = str::length(string);
  const int offset = strlen - str::length(suffix);
  if (offset < 0)
    return false;
  return std::equal(str::begin(string) + offset, str::begin(string) + strlen, str::begin(suffix));
}

// Modifies the specified C-style string such that it begins at beginIdx and ends at endIdx - 1.
inline void substr(char*& string, const int beginIdx, const int endIdx) {
  assert(0 <= beginIdx); assert(beginIdx <= endIdx); assert(endIdx <= std::strlen(string));
  string[endIdx] = '\0';
  string += beginIdx;
}
