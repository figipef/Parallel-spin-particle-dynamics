#pragma once

#include <vector>
#include <iostream>

// Custom struct to hold either a 1D vector or a 2D matrix
struct Histogram {
    bool is_matrix;  // Flag to determine if it's a matrix
    std::vector<int> vector1D;
    std::vector<std::vector<int>> matrix2D;

    // Flatten the histogram for MPI communication
    std::vector<int> flatten() const {
        if (!is_matrix) {
            return vector1D;
        } else {
            std::vector<int> flat;
            for (const auto& row : matrix2D) {
                flat.insert(flat.end(), row.begin(), row.end());
            }
            return flat;
        }
    }

    // Resize based on bin sizes
    void resize(int* bin_n) {
        if (!is_matrix) {
            vector1D.resize(bin_n[0], 0);
        } else {
            matrix2D.resize(bin_n[1], std::vector<int>(bin_n[0], 0));
        }
    }

    // Restore from flattened data
    void from_flattened(const std::vector<int>& flat, int* bin_n) {
        if (!is_matrix) {
            vector1D = flat;
        } else {
            int idx = 0;
            for (int i = 0; i < bin_n[1]; ++i) {
                for (int j = 0; j < bin_n[0]; ++j) {
                    matrix2D[i][j] = flat[idx++];
                }
            }
        }
    }
};

Histogram createHistogram(const int, const int*);
