#pragma once

#include "particle.hpp"
#include <vector>
#include <fstream> // Removable
#include <iostream>
#include <cmath> 

int test(int);

double* cross(double[3], double[3]);
double inner(double[3], double[3]);
double gamma(double[3]);
void boris(Particle*, double, double[3], double[3], int);
void createParticles(Particle*, int);
void PerformDiagnostics(std::vector<int>*&, Particle*,  \
	int , std::string , std::string , std::string , double , double , double, \
	double , double , double , double , double , double );
void PerformDiagnostics(std::vector<int>*&, Particle,std::string*, double*, double*, double*, int);