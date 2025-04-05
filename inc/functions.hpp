#pragma once

#include "particle.hpp"
#include "laser.hpp"
#include "histogram.hpp"
#include "diagnostics.hpp"
#include "writers.hpp"

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <random>

double* cross(double[3], double[3]); // Calculates the cross product between 2 vectors
double* Omega(double, double[3], double[3], double[3]); // Calculates the Omega used in spin evolution
double inner(double[3], double[3]); // Calculates the inner product between 2 vectors
double gamma(double[3]); // Calculates the gamma for a vector
double* f_RR(double, double[3], double[3], double[3]); // Calculates the radiation reaction
void PerformDiagnostics(Histogram&, Particle,std::string*, double*, double*, double*, int); // Writes the diagnostics to the files
// Particle pusher (BORIS + RR + BMT)
void boris(Particle* particles, Laser* lasers, double time, double time_step, int n_of_particles, int n_of_lasers, int RR, Histogram* hist = nullptr, DiagnosticParameters* diag_params = nullptr);
// Creates the particles
void createParticles(Particle*, int, std::string*, double*, double*, double*, double*, Laser*, int);
// Writes the fields to perform diagnsotics
void FieldDiagWritter(double&, int&, double*& fieldiag, Laser*& lasers);
// Sets up the input variables
void setupInputVariable(std::ifstream&, int&, std::string*&, double*&, double*&, double*&, double*&, double&, double&, int&, std::string*&, double*&, double*&, double*&, int*&, int&, double*&, Laser*&, int&, int&);
// Sets up the parameters to follow and save the particle data to files
void setupFollowParticles(int*, std::ofstream*, std::ofstream*, std::ofstream*, int*&, int&, int);