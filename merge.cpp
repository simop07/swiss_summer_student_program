// Compile using "g++ -fsanitize=address merge.cpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

constexpr int nMinAnalysedRows{6};  // Minimum index of analysed rows INCLUDED
constexpr int nMaxAnalysedRows{8008};  // Maximum rows EXCLUDED

void mergeCsv() {
  // Create file to write on
  std::ofstream outputFile{"./data/4LayersD.txt"};

  // Read input files
  std::string fileRoot{"./data/4LayerD/Autosave/C1SCOPE#"};
  std::string fileType{".csv"};
  std::string line{};

  // Loop on files to merge
  for (int i = 200001; i <= 201011; ++i) {
    std::cout << "\n*** FILE n. " << i << " ***\n";
    std::string fileName{fileRoot + std::to_string(i) + fileType};
    std::ifstream currentInfile(fileName);
    int row{1};

    // Check if file exists
    if (!currentInfile.is_open()) {
      std::cerr << "Error: could not be able to open file " << fileName << '\n';
      return;
    }

    // Loop on rows for each file
    while (std::getline(currentInfile, line)) {
      if (row < nMinAnalysedRows) {
        ++row;
        continue;
      }
      if (row >= nMaxAnalysedRows) {
        break;
      }

      // Defining loop variables
      std::stringstream ss(line);
      std::string item{};
      int column{1};

      // std::cout << "\n*** ROW n. " << row << " ***\n";

      // Loop on columns for each file
      while (std::getline(ss, item, ',')) {
        if (item.empty()) {
          continue;
        }

        if (column == 1 && row == 6) {
          // std::cout << "Timestamp = " << std::stod(item) << '\n';
          outputFile << "0\t0\t" << std::stod(item) << "\t0\t0\t0\t";
        } else if (column == 2) {
          // std::cout << "Amplitude = " << std::stod(item) << '\n';
          outputFile << std::stod(item) << '\t';
        }
        ++column;
      }
      ++row;
    }
    outputFile << '\n';
  }
  outputFile.close();
}

int main() {
  mergeCsv();

  return EXIT_SUCCESS;
}