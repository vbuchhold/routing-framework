#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

#include <csv.h>

#include "DataStructures/Geometry/CoordinateTransformation.h"
#include "Tools/CommandLine/CommandLineParser.h"
#include "Tools/Math.h"
#include "Tools/StringHelpers.h"

inline void printUsage() {
  std::cout <<
      "Usage: VisumExtractPoi -crs <file> -i <file> -o <file>\n"
      "Extracts points of interest (e.g. schools and hospitals) from given Visum files.\n"
      "  -crs <file>       coordinate reference system used in the Visum network files\n"
      "  -i <file>         directory that contains the Visum network files\n"
      "  -o <file>         place output in <file>\n"
      "  -help             display this help and exit\n";
}

int main(int argc, char* argv[]) {
  try {
    CommandLineParser clp(argc, argv);
    if (clp.isSet("help")) {
      printUsage();
      return EXIT_SUCCESS;
    }

    const auto crs = clp.getValue<int>("crs");
    const auto visumPathName = clp.getValue<std::string>("i");
    auto outputFileName = clp.getValue<std::string>("o");
    if (!endsWith(outputFileName, ".csv"))
      outputFileName += ".csv";

    CoordinateTransformation trans(crs, CoordinateTransformation::WGS_84);
    std::ofstream outputFile(outputFileName);
    if (!outputFile.good())
      throw std::invalid_argument("file cannot be opened -- '" + outputFileName + "'");
    outputFile << std::fixed;
    outputFile << "# Source: " << visumPathName << "\n";
    outputFile << "category,longitude,latitude\n";

    std::cout << "Reading POI categories..." << std::flush;
    std::map<int, std::string> poiCategories;
    int id;
    char* name;
    using VisumFileReader = io::CSVReader<2, io::trim_chars<>, io::no_quote_escape<';'>>;
    VisumFileReader poiCatFileReader(visumPathName + "/POIKATEGORIE.csv");
    poiCatFileReader.read_header(io::ignore_extra_column, "NR", "NAME");
    while (poiCatFileReader.read_row(id, name)) {
      assert(*name != '\0');
      assert(poiCategories.find(id) == poiCategories.end());
      poiCategories[id] = name;
    }
    std::cout << " done.\n";

    for (const auto& cat : poiCategories) {
      const auto poiFileName = visumPathName + "/POIOFCAT_" + std::to_string(cat.first) + ".csv";
      std::ifstream poiFile(poiFileName);
      if (poiFile.good()) {
        std::cout << "Extracting POIs of category " << cat.second << "..." << std::flush;
        double easting, northing;
        double lng, lat;
        VisumFileReader poiFileReader(poiFileName, poiFile);
        poiFileReader.read_header(io::ignore_extra_column, "XKOORD", "YKOORD");
        while (poiFileReader.read_row(easting, northing)) {
          trans.forward(easting, northing, lng, lat);
          outputFile << cat.second << ',' << toDegrees(lng) << ',' << toDegrees(lat) << '\n';
        }
        std::cout << " done.\n";
      }
    }
  } catch (std::exception& e) {
    std::cerr << argv[0] << ": " << e.what() << '\n';
    std::cerr << "Try '" << argv[0] <<" -help' for more information.\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
