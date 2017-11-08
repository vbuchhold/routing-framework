#pragma once

#include <cassert>
#include <stdexcept>
#include <string>
#include <unordered_map>

// A facility for translating strings into enum values. Enums must specialize initNameToEnumMap.
template <typename T>
class EnumParser {
 public:
  // Constructs and initializes an enum parser.
  EnumParser() {
    initNameToEnumMap();
  }

  // Fills the map with name/value pairs. Enumerations must specialize this member function.
  void initNameToEnumMap() {
    assert(false);
  }

  // Returns the enum value with the specified name.
  T operator()(const std::string& name) const {
    const auto iter = nameToEnum.find(name);
    if (iter == nameToEnum.end())
      throw std::invalid_argument("no enum value with the specified name -- '" + name + "'");
    return iter->second;
  }

 private:
  std::unordered_map<std::string, T> nameToEnum; // A map to translate strings into enum values.
};
