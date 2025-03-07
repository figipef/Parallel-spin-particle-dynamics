#pragma once

#include "particle.hpp"
#include "laser.hpp"
#include <vector>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <sstream>
#include <cmath> 

double* cross(double[3], double[3]); // Calculates the cross product between 2 vectors
double* Omega(double, double[3], double[3], double[3]); // Calculates the Omega used in spin evolution
double inner(double[3], double[3]); // Calculates the inner product between 2 vectors
double gamma(double[3]); // Calculates the gamma for a vector
void boris(Particle*&, Laser*, double, double,  int, int); // Simple boris pusher
void createParticles(Particle*, int); // Creates particles
void PerformDiagnostics(std::vector<int>*&, Particle,std::string*, double*, double*, double*, int);
void setupInputVariable(std::ifstream&, int&, double&, double&, std::string*&, double*&, double*&, double*&, double*&, int&, Laser*&, int&);
void writeToFile(std::ofstream&,  const Particle&, char);
void parseVector(const std::string&, double[3]);