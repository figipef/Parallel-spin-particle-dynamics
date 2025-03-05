#pragma once

#include <cmath> 
#include <algorithm>

class Laser
{
private:

	double e_0[3]; // Eletric field vector
	double b_0[3]; // Magnetic field vector
	double k[3]; // Wave vector
	double freq; // frequency calculated from wave vector
	double env_freq; // frequency for the envelope
	double length; // length of the envelope
	int tag; // Tag for identification
	int env;

public:

	// default constructor
	Laser();
	// constructor
	Laser(double[3], double[3], int); // intensity vector (Eletric), wavevector, tag
	// constructor for an envelope
	Laser(double[3], double[3], int, double, double); // intensity vector (Eletric), wavevector, tag
	
	double* get_E(double[3], double); // Return the true electric field at a point at a time
	double* get_B(double[3], double); // Return the true magnetic field at a point at a time

	double* get_time_E(double); // Return the electric field independent of positon
	double* get_time_B(double); // Return the magnetic field independent of positon

	const double* get_E_0() const; // Return the original electric field (CONSTANT FIELD)
	const double* get_B_0() const; // Return the original magnetic field (CONSTANT FIELD)

	void set_B_0(double[3]);

};