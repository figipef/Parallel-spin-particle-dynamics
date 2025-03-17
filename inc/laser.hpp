#pragma once

#define _USE_MATH_DEFINES

#include <cmath> 
#include <algorithm>
#include <iostream>

class Laser
{
private:

	double e_0[3]; // Eletric field vector
	double b_0[3]; // Magnetic field vector
public:
	double k[3]; // Wave vector
	double freq; // frequency calculated from wave vector
	double env_freq; // frequency for the envelope
	double length; // length of the envelope
	int tag; // Tag for identification
	int type;
	int ext_phase;

	// default constructor
	Laser();
	// constructor with basica structure for an eletromagnetic field
	Laser(double[3], double[3], int, int _ext_phase = 0); // intensity vector (Eletric), wavevector, tag
	// constructor with the option of b or k field and on/off frequency
	Laser(double[3], double[3], int, char, char, int _ext_phase = 0, double* _freq = nullptr); // intensity vector (Eletric), wavevector, tag (Para usar frequência, pôr &_freq)
	// constructor for an envelope
	Laser(double[3], double[3], int, double, double, int _ext_phase = 0); // intensity vector (Eletric), wavevector, tag and envelope parameters
	
	double* get_E(double[3], double); // Return the true electric field at a point at a time
	double* get_B(double[3], double); // Return the true magnetic field at a point at a time

	double* get_time_E(double); // Return the electric field independent of positon
	double* get_time_B(double); // Return the magnetic field independent of positon

	const double* get_E_0() const; // Return the original electric field (CONSTANT FIELD)
	const double* get_B_0() const; // Return the original magnetic field (CONSTANT FIELD)

	void set_B_0(double[3]);

	double* get_fields_envelope(double[3], double); // Return an array with the Electric and Magnetic Fields (E,B) for the envelope defined

	// GENERAL GET FIELDS FUNCTION!!!!!!! Should be used to get the fields of a given laser
	double* get_fields(double* t = nullptr, const double pos[3] = nullptr, int* _type = nullptr); // Para usar o tempo, pôr &t!!!!
};