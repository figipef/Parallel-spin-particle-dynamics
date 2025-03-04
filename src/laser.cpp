#include "laser.hpp"

Laser::Laser(){} // Default Constructor

Laser::Laser(double _e_0[3], double _k[3], int _tag) : tag(_tag){
	std::copy(_e_0, _e_0 + 3, e_0); // copy the letric field intensity from _e_0 to e_0
	std::copy(_k, _k+3, k); // copy the wave vector _k to k

	b_0[0] = (k[1]*e_0[2] - k[2]*e_0[1])/1.; // calculate the cross product the  must be divided by c=1
	b_0[1] = (k[2]*e_0[0] - k[0]*e_0[2])/1.;
	b_0[2] = (k[0]*e_0[1] - k[1]*e_0[0])/1.;

	freq = 1. * sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]); // calculate the frequency based on the wave vector

}

double* Laser::get_E(double pos[3], double t){

	double* e_field = new double[3];
	double aux = cos(pos[0] * k[0] + pos[1] * k[1] + pos[2] * k[2]- freq * t);

	for (int i = 0; i < 3; i++){
		e_field[i] = e_0[i] * aux;
	}

	return e_field;
} 
double* Laser::get_B(double pos[3], double t){

	double* b_field = new double[3];
	double aux = cos(pos[0] * k[0] + pos[1] * k[1] + pos[2] * k[2]- freq * t);
	
	for (int i = 0; i < 3; i++){
		b_field[i] = b_0[i] *aux;
	}

	return b_field;

} 

double* Laser::get_time_E(double t){

	double* e_field = new double[3];

	for (int i = 0; i < 3; i++){
		e_field[i] = e_0[i] * cos(freq * t);
	}

	return e_field;

} 

double* Laser::get_time_B(double t){

	double* b_field = new double[3];

	for (int i = 0; i < 3; i++){
		b_field[i] = b_0[i] * cos(freq * t);
	}

	return b_field;

} 

double* Laser::get_E_0(){
	return e_0;
} 

double* Laser::get_B_0(){
	return b_0;
} 

void Laser::set_B_0(double new_B[3]){
	std::copy(new_B, new_B + 3, b_0);
}