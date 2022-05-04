#ifndef COMMONMACROS_H
#define COMMONMACROS_H

#include <vector>
#include <string.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <omp.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "Eigen/Eigen"

#include <chrono>
#include <thread>

#define LOOP_l_N(NUM)  for(int l = 0; l < NUM ; ++l)
#define LOOP_k_N(NUM)  for(int k = 0; k < NUM ; ++k)

static std::vector<double> Tokenize(const std::string stringIn,
                                               const char *delimiters) {
  // Member function ID string
  std::string memberID = "Tokenize";

  std::vector<double> outputVector;

  // flowControl while (for error handling)
  int flowControl = 1;
  while (flowControl == 1) {
    // Tokenize string read from file (spaces or commas are allowed(
    char *token = strtok((char *)stringIn.c_str(), delimiters);
    while (token != nullptr) {
      outputVector.push_back(std::atof(token));
      token = strtok(nullptr, delimiters);
    }

    // Continue flow normally (no unrecoverable errors up to this point)
    flowControl++;
  }

  return outputVector;
}

#endif // COMMONMACROS_H
