#pragma once

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "Tools/Constants.h"

// An origin-destination (OD) pair, representing a travel demand or a query.
struct OriginDestination {
  // Constructs an OD-pair from o to d.
  OriginDestination(const int o, const int d)
      : origin(o), destination(d) {}

  // Compares this OD-pair with rhs lexicographically.
  bool operator<(const OriginDestination& rhs) const {
    return origin < rhs.origin || (origin == rhs.origin && destination < rhs.destination);
  }

  int origin;
  int destination;
};

// An origin-destination (OD) pair that additionally stores an origin zone and a destination zone.
// Zones or traffic cells represent for example residential or commercial areas.
struct ClusteredOriginDestination : public OriginDestination {
  // Constructs a clustered OD-pair from o to d.
  ClusteredOriginDestination(const int o, const int d, const int oZone, const int dZone)
      : OriginDestination(o, d), originZone(oZone), destinationZone(dZone) {}

  // Compares this clustered OD-pair with rhs lexicographically.
  bool operator<(const ClusteredOriginDestination& rhs) const {
    if (originZone < rhs.originZone)
      return true;
    if (originZone > rhs.originZone)
      return false;
    if (destinationZone < rhs.destinationZone)
      return true;
    if (destinationZone > rhs.destinationZone)
      return false;
    return OriginDestination::operator<(rhs);
  }

  int originZone;
  int destinationZone;
};

// Reads the specified file into a vector of OD-pairs.
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
        line >> lineType >> o >> d;
        assert(line);
        assert(o >= 0);
        assert(d >= 0);
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

// Reads the specified file into a vector of clustered OD-pairs.
std::vector<ClusteredOriginDestination> importClusteredODPairsFrom(std::ifstream& in) {
  std::vector<ClusteredOriginDestination> pairs;
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
        int o, d, oZone = INVALID_ID, dZone = INVALID_ID;
        line >> lineType >> o >> d >> std::ws;
        if (!line.eof()) {
          line >> oZone >> dZone >> std::ws;
          assert(oZone >= 0);
          assert(dZone >= 0);
        }
        assert(o >= 0);
        assert(d >= 0);
        assert(line); assert(line.eof());
        pairs.emplace_back(o, d, oZone, dZone);
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
