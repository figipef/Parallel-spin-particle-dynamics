#pragma once

// Struct to help passing the diagnostic information to functions
struct DiagnosticParameters {
    std::string* params;     // Pointer to array of strings
    double* bsize;           // Pointer to array of doubles (size values)
    double* bmax;            // Pointer to array of doubles (max values)
    double* bmin;            // Pointer to array of doubles (min values)
    int n_of_pars;           // Number of parameters (int)
    
    // Constructor for easy initialization
    DiagnosticParameters(std::string* p, double* bs, double* bm, double* bn, int n)
        : params(p), bsize(bs), bmax(bm), bmin(bn), n_of_pars(n) {}
};
