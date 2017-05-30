#pragma once

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// An origin-destination (OD) pair, representing a travel demand or a query.
struct OriginDestination {
  // Constructs an OD-pair from origin to destination.
  OriginDestination(const int origin, const int destination)
      : origin(origin), destination(destination) {}

  int origin;
  int destination;
};

// Reads the specified file containing OD-pairs into a vector.
std::vector<OriginDestination> importODPairsFrom(std::ifstream& in) {
  std::vector<OriginDestination> pairs;
  char lineType;
  std::string lineStr;

  // Read the file line by line.
  while (getline(in, lineStr)) {
    assert(!lineStr.empty());
    switch (lineStr[0]) {
      case 'c':
      case 'p':
      case 'r':
      case 'd': {
        // Skip all non-OD-pair lines (e.g., comment lines).
        break;
      }
      case 'q': {
        // Read an OD-pair line.
        std::istringstream line(lineStr);
        int o, d;
        line >> lineType >> o >> d >> std::ws;
        assert(o >= 0);
        assert(d >= 0);
        assert(line); assert(line.eof());
        pairs.emplace_back(o, d);
        break;
      }
      default: {
        assert(false);
      }
    }
  }

  assert(in.eof());
  return pairs;
}
