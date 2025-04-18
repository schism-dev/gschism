#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <string>
#include <stdexcept>
// #include <format>
#include "ansi.hpp"

namespace py = pybind11;

const double pi = 3.14159265358979323846;
const double deg2rad = pi / 180.0;

// Helper function for computing the local node index at a given offset from another.
inline int nxq(int offset, int index, int size) {
  assert(index + offset + 10*size >= 0);
  // In c++, the % operator can return negative numbers, for this reason we add 10*size 
  // to the index to make it positive
  return (index + offset + 10*size) % size; 
}

// ----------------------------------------------------
// IO functions

int readInt(std::ifstream& file) {
  // Read an integer from the file. The rest of the line is ignored.

  std::string line;
  std::getline(file, line);
  try {
    std::istringstream iss(line);
    int value;
    iss >> value;
    return value;

  } catch (const std::exception& e) {
    std::cerr << "Error reading integer from line: " << line << std::endl;
    throw;
  }
}

std::vector<int> readInts(std::ifstream& file) {
  // Read a list of integers from the file. The rest of the line is ignored.
  std::string line;
  std::getline(file, line);
  printf("line = %s\n", line.c_str());
  try {
    std::istringstream iss(line);
    std::vector<int> values;
    int value;
    while (iss >> value) {
      values.push_back(value);
    }
    return values;
  } catch (const std::exception& e) {
    std::cerr << "Error reading integer from line: " << line << std::endl;
    throw;
  }
}

// ----------------------------------------------------
// Printing, logging, etc.

void print(const std::string& message) {
  // Print a message to the console.

  py::print(message);
}

void printFormatted(const std::string& message, ...) {
  // Print a formatted message to the console.

  char buffer[2048];
  va_list args;
  va_start(args, message);
  vsnprintf(buffer, sizeof(buffer), message.c_str(), args);
  va_end(args);
  py::print(std::string(buffer));
}

void errorMessage(const std::string& message, ...) {
  // Print an error message to the console.

  char buffer[2048];
  va_list args;
  va_start(args, message);
  vsnprintf(buffer, sizeof(buffer), message.c_str(), args);
  va_end(args);
  py::print(std::string(buffer));
}

void progressBar(const std::string& desc, int current, int total, const std::string& color, int barWidth = 20) {
  // Print a progress bar to the console, tqdm style. Colors are found in ansi.hpp.
  // The bar is updated only if the percentage has changed.
    static int currentPercent = 0;

    total = (total <= 1)?2:total; // Edge case where total is 1.

    int percent = static_cast<int>(current * 100.0 / (total - 1));
    if (percent != currentPercent) {
        py::print(ansi::progressBar(desc, current, total - 1, color, barWidth));
        py::print(ansi::cursorUp(2));
    }
    currentPercent = percent;
}
