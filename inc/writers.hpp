#pragma once

#include "histogram.hpp"
#include "particle.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

void writeToFile(std::ofstream&,  const Particle&, char); // To follow a single particle and its trajectory
void writeDiagnosticsToFile(const Histogram, const int, const double, const std::string*, const int); // To write the Diagnostics to a numbered file
void parseVector(const std::string&, double[3]); // To make a vector a string
void printProgressBar(int, int, int);