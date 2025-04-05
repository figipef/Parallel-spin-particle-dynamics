#include "laser.hpp"

Laser::Laser(){} // Default Constructor

Laser::Laser(double _e_0[3], double _k[3], int _tag, int _ext_phase) : tag(_tag), ext_phase(_ext_phase){
	std::copy(_e_0, _e_0 + 3, e_0); // copy the letric field intensity from _e_0 to e_0
	std::copy(_k, _k+3, k); // copy the wave vector _k to k

	b_0[0] = (k[1]*e_0[2] - k[2]*e_0[1])/1.; // calculate the cross product the  must be divided by c=1
	b_0[1] = (k[2]*e_0[0] - k[0]*e_0[2])/1.;
	b_0[2] = (k[0]*e_0[1] - k[1]*e_0[0])/1.;

	freq = 1. * sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]); // calculate the frequency based on the wave vector

	type = 2; // type is defined to know if it is initalized
}

Laser::Laser(double _e_0[3], double _temp[3], int _tag, char _f_option, char _t_option, int _ext_phase, double* _freq) : tag(_tag), ext_phase(_ext_phase){
	// constructor for a bunch of stuff
	std::copy(_e_0, _e_0 + 3, e_0); // copy the letric field intensity from _e_0 to e_0

	int no_k = 1;

	if (_f_option == 'k' || _f_option == 'K'){

		std::copy(_temp, _temp+3, k); // copy the wave vector _k to k

		b_0[0] = (k[1]*e_0[2] - k[2]*e_0[1])/1.; // calculate the cross product the  must be divided by c=1
		b_0[1] = (k[2]*e_0[0] - k[0]*e_0[2])/1.;
		b_0[2] = (k[0]*e_0[1] - k[1]*e_0[0])/1.;

		no_k = 0;

	} else if (_f_option == 'b' || _f_option == 'B'){

		std::copy(_temp, _temp+3, b_0);

		k[0] = 0;
		k[1] = 0;
		k[2] = 0;

	} else {throw std::runtime_error("Not a valid field option");}

	if (_t_option == 'y' || _t_option == 'Y'){
		type = 1;
		if (no_k || _freq){
			freq = *_freq;
		} else if (!no_k){
			freq = 1. * sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]); // calculate the frequency based on the wave vector c = 1
		} else {throw std::runtime_error("No k-vector or frequency given to calculate frequency");}

	}else if (_t_option == 'n' || _t_option == 'N'){

		freq = 0;
		type = 0;

	}else {throw std::runtime_error("Not a valid time option");}
}

Laser::Laser(double _e_0[3], double _k[3], int _tag, double _length, int _ext_phase) : length(_length), tag(_tag), ext_phase(_ext_phase){
	std::copy(_e_0, _e_0 + 3, e_0); // copy the letric field intensity from _e_0 to e_0
	
	//std::cout <<"E" << e_0[0]<<e_0[1]<<e_0[2]<< "\n ";
	std::copy(_k, _k+3, k); // copy the wave vector _k to k

	b_0[0] = (k[1]*e_0[2] - k[2]*e_0[1])/1.; // calculate the cross product the  must be divided by c=1
	b_0[1] = (k[2]*e_0[0] - k[0]*e_0[2])/1.;
	b_0[2] = (k[0]*e_0[1] - k[1]*e_0[0])/1.;

	freq = 1. * sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]); // calculate the frequency based on the wave vector c = 1

	type = 3;
}

const double* Laser::get_E_0() const {
	return e_0;
} 

const double* Laser::get_B_0() const {
	return b_0;
} 

void Laser::set_B_0(double new_B[3]){
	std::copy(new_B, new_B + 3, b_0);
}

double* Laser::get_fields(double* t, const double pos[3], int* _type){

	if (_type){
		type = *_type;
	}

	if (type == 0){

		// Return the constant fields
		double* fields = new double[6];
		for (int i = 0; i < 3; ++i) fields[i] = e_0[i];
    	for (int i = 3; i < 6; ++i) fields[i] = b_0[i - 3];
    	return fields;

	} else if (type == 1){

		// Return oscilating fields equal in all space
		if (!t){throw std::runtime_error("No time given, but type = 1");}
		if (freq == 0){throw std::runtime_error("No frequency given, but type = 1");}

		double* fields = new double[6];

		double aux = cos(freq * (*t) + ext_phase * M_PI / 2);

		for (int i = 0; i < 3; ++i) fields[i] = e_0[i] * aux;
    	for (int i = 3; i < 6; ++i) fields[i] = b_0[i - 3] * aux;
    	return fields; 

	} else if (type == 2){

		// Return the laser fields
		if (!t){throw std::runtime_error("No time given, but type = 2");}
		if (!pos){throw std::runtime_error("No position given, but type = 2");}
		if (freq == 0){throw std::runtime_error("No frequency given, but type = 2");}

		double* fields = new double[6]; // E[3] x B[3] an array with the eletric and magnetic fields

		double xdk = pos[0] * k[0] + pos[1] * k[1] + pos[2] * k[2]; // k.x for future calculations

		double field_osc = cos(xdk - freq * (*t) + ext_phase * M_PI / 2);// The field part responsible for oscilation

		for (int i = 0; i < 3; i++) fields[i] = e_0[i] * field_osc;
		for (int i = 3; i < 6; i++) fields[i] = b_0[i - 3] * field_osc;
		

		return fields;

	} else if (type == 3){

		if (!t){throw std::runtime_error("No time given, but type = 3");}
		if (!pos){throw std::runtime_error("No position given, but type = 3");}
		if (freq == 0){throw std::runtime_error("No frequency given, but type = 3");}

		double* fields = new double[6]; // E[3] x B[3] an array with the eletric and magnetic fields

		double xdk = pos[0] * k[0] + pos[1] * k[1] + pos[2] * k[2]; // k.x for future calculations
		//std::cout <<"xdk = "<< (xdk - freq * (*t)) / (length * M_PI)<<"\n";
		double field_osc = cos(xdk - freq * (*t) + ext_phase * M_PI / 2);// The field part responsible for oscilation

 		// Calculate the phase normalized to the length so that the function is 0 at phase = 0 and phase = L
		double phase = (xdk - freq * (*t))* M_PI / length ;
		//std::cout <<" phase: "<< phase <<"\n";
		double envelope = 0;

		if (phase > 0 && phase < M_PI){
			// Change the sin * sin to other functions or to a different function that has to be defined by a character probably
			envelope = sin(phase) * sin(phase); // Calculate the value of the phase | Envelope function
		} else {
			envelope = 0;
		}

		for (int i = 0; i < 3; i++) fields[i] = e_0[i] * field_osc * envelope;
		
		for (int i = 3; i < 6; i++) fields[i] = b_0[i - 3] * field_osc * envelope;

		return fields;

	} else {throw std::runtime_error("Not a valid type");}
}