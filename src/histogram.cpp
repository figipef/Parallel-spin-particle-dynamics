#include "histogram.hpp"

Histogram createHistogram(const int params, const int* bin_n) {
    
    // Function to create an histogram according to its dimensions

    if (params > 2){
        throw std::runtime_error("createHistogram: Number of diagnostic parameters larger than 2");
    }

    Histogram hist;

    if (params == 1) { // Create a single vector inside
        hist.is_matrix = false;
        hist.vector1D = std::vector<int>(bin_n[0], 0);
    } 
    else if (params == 2) { // Create a matrix according to the dimensions of the parameters
        hist.is_matrix = true;
        hist.matrix2D = std::vector<std::vector<int>>(bin_n[1], std::vector<int>(bin_n[0], 0));
    }

    return hist;
}
