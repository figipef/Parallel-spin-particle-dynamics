#pragma once

#define _USE_MATH_DEFINES

#include <cmath> 
#include <algorithm>
#include <iostream>

class Laser
{

public:
	
	double e_0[3]; // Eletric field vector
	double b_0[3]; // Magnetic field vector
	double k[3]; // Wave vector
	double freq; // frequency calculated from wave vector
	double length; // length of the envelope
	int tag; // Tag for identification
	int type; // Type of Laser
	int ext_phase; // external phase (increments of PI/2)

	// default constructor
	Laser();
	// constructor with basica structure for an eletromagnetic field
	Laser(double[3], double[3], int, int _ext_phase = 0); // intensity vector (Eletric), wavevector, tag
	// constructor with the option of b or k field and on/off frequency
	Laser(double[3], double[3], int, char, char, int _ext_phase = 0, double* _freq = nullptr); // intensity vector (Eletric), wavevector, tag (Para usar frequência, pôr &_freq)
	// constructor for an envelope
	Laser(double[3], double[3], int, double, int _ext_phase = 0); // intensity vector (Eletric), wavevector, tag and envelope parameters

	const double* get_E_0() const; // Return the original electric field (CONSTANT FIELD)
	const double* get_B_0() const; // Return the original magnetic field (CONSTANT FIELD)

	void set_B_0(double[3]);

	// GENERAL GET FIELDS FUNCTION!!!!!!! Should be used to get the fields of a given laser
	double* get_fields(double* t = nullptr, const double pos[3] = nullptr, int* _type = nullptr); // Para usar o tempo, pôr &t!!!!
};