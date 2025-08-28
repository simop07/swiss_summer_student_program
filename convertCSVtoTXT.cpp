#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

void convertCSVtoTXT(std::string const& folderPath) {
  // Loop over files of directory
  for (auto const& fileCSV : std::filesystem::directory_iterator(folderPath)) {
    // Select file validity and .CSV extension
    if (fileCSV.is_regular_file() && fileCSV.path().extension() == ".CSV") {
      // Input file name opening
      std::ifstream infile(fileCSV.path());
      if (!infile.is_open()) {
        std::cerr << "Failed to open " << fileCSV.path() << '\n';
        return;
      }

      // Replace extension from .csv to .txt
      std::filesystem::path txtPath = fileCSV.path();
      txtPath.replace_extension(".txt");

      // Skip to next file if already converted
      if (std::filesystem::exists(txtPath)) {
        std::cout << "Skipping (already converted): " << fileCSV.path() << '\n';
        continue;
      }

      // Create output .txt file
      std::ofstream outFile(txtPath);
      if (!outFile.is_open()) {
        std::cerr << "Failed to create " << txtPath << '\n';
        return;
      }

      // Read input file lines
      std::string line;
      while (std::getline(infile, line)) {
        // Defining loop variables
        std::stringstream ss(line);
        std::string item;
        bool first = true;

        // Loop over columns
        while (std::getline(ss, item, ';')) {
          if (!first) {
            outFile << '\t';
          }
          outFile << item;
          first = false;
        }
        outFile << '\n';
      }
      outFile.close();

      // Status
      std::cout << "Converted: " << fileCSV.path() << " -> " << txtPath << '\n';
    }
  }
}

int main() {
  convertCSVtoTXT("/mnt/c/Users/Simone/Desktop/ConversionFolder");

  return EXIT_SUCCESS;
}