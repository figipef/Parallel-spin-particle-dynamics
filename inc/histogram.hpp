#pragma once

#include <vector>
#include <iostream>

// Custom struct to hold either a 1D vector or a 2D matrix
struct Histogram {
    bool is_matrix;  // Flag to determine if it's a matrix
    std::vector<int> vector1D;
    std::vector<std::vector<int>> matrix2D;
};

Histogram createHistogram(const int, const int*);
