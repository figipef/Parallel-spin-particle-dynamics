#pragma once

#include "particle.hpp"
#include "laser.hpp"
#include <vector>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <sstream>
#include <cmath> 

int test(int);

double* cross(double[3], double[3]);
double* Omega(double, double[3], double[3], double[3]);
double inner(double[3], double[3]);
double gamma(double[3]);
void boris(Particle*&, Laser*, double,  int, int);
void createParticles(Particle*, int);
void PerformDiagnostics(std::vector<int>*&, Particle*,  \
	int , std::string , std::string , std::string , double , double , double, \
	double , double , double , double , double , double );
void PerformDiagnostics(std::vector<int>*&, Particle,std::string*, double*, double*, double*, int);
void setupInputVariable(std::ifstream&, int&, double&, double&, std::string*&, double*&, double*&, double*&, double*&, int&);
void writeToFile(std::ofstream&,  const Particle&, char);