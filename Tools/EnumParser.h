#pragma once

#include <stdexcept>
#include <string>
#include <unordered_map>

// A facility for translating strings into enum values. Enumerations should specialize the default
// constructor and fill the map in the constructor.
template <typename T>
class EnumParser {
 public:
  // Default constructor to be specialized.
  EnumParser() {}

  // Returns the enum value with the specified name.
  T parse(const std::string& name) const {
    const auto iter = nameToEnum.find(name);
    if (iter == nameToEnum.end())
      throw std::invalid_argument("no enum value with the specified name -- '" + name + "'");
    return iter->second;
  }

 private:
  std::unordered_map<std::string, T> nameToEnum; // A map to translate strings into enum values.
};
