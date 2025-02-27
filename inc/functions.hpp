#pragma once

#include "particle.hpp"
#include <fstream>

int test(int);

double* cross(double[3], double[3]);
double gamma(double[3]);
void boris(Particle*, double, double[3], double[3], int);
void createParticles(Particle*, int);
void PerformDiagnostics(std::ifstream&, Particle*,  \
	int , std::string , std::string , std::string , double , double , double, \
	double , double , double , double , double , double );