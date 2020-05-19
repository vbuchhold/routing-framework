#pragma once

#include <cassert>
#include <string>
#include <unordered_map>

#include "Tools/Logging/NullLogger.h"

// This class manages a set of named loggers. The logger type should be std::basic_ofstream or
// NullLogger. In the former case, each logger writes to a file whose name is the concatenation of a
// common base file name with the name of the logger. In the latter case, each logger discards the
// data written to it, avoiding any overhead at runtime. When a logger is requested, we return a
// reference to it if it already exists. Otherwise, the logger is newly constructed, an optional
// header line is written to it, and a reference to it is returned.
template <typename LoggerT = NullLogger>
class LogManager {
 public:
  LogManager() = delete;

  static LoggerT& getLogger(const std::string& name, const std::string& header = "") {
    auto iter = loggers.find(name);
    if (iter == loggers.end()) {
      iter = loggers.emplace(name, LoggerT(baseFileName + name)).first;
      assert(iter->second);
      iter->second << header;
    }
    return iter->second;
  }

  static void setBaseFileName(const std::string& fileName) {
    baseFileName = fileName;
  }

 private:
  inline static std::unordered_map<std::string, LoggerT> loggers;
  inline static std::string baseFileName;
};
