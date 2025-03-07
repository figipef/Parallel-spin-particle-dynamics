#include "laser.hpp"

Laser::Laser(){} // Default Constructor

Laser::Laser(double _e_0[3], double _k[3], int _tag) : tag(_tag){
	std::copy(_e_0, _e_0 + 3, e_0); // copy the letric field intensity from _e_0 to e_0
	std::copy(_k, _k+3, k); // copy the wave vector _k to k

	b_0[0] = (k[1]*e_0[2] - k[2]*e_0[1])/1.; // calculate the cross product the  must be divided by c=1
	b_0[1] = (k[2]*e_0[0] - k[0]*e_0[2])/1.;
	b_0[2] = (k[0]*e_0[1] - k[1]*e_0[0])/1.;

	freq = 1. * sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]); // calculate the frequency based on the wave vector

	type = 2; // type is defined to know if it is initalized
}

Laser::Laser(double _e_0[3], double _temp[3], int _tag, char _f_option, char _t_option, double* _freq) : tag(_tag){
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

Laser::Laser(double _e_0[3], double _k[3], int _tag, double _length, double _env_freq) : tag(_tag), length(_length), env_freq(_env_freq){
	std::copy(_e_0, _e_0 + 3, e_0); // copy the letric field intensity from _e_0 to e_0
	std::copy(_k, _k+3, k); // copy the wave vector _k to k

	b_0[0] = (k[1]*e_0[2] - k[2]*e_0[1])/1.; // calculate the cross product the  must be divided by c=1
	b_0[1] = (k[2]*e_0[0] - k[0]*e_0[2])/1.;
	b_0[2] = (k[0]*e_0[1] - k[1]*e_0[0])/1.;

	freq = 1. * sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]); // calculate the frequency based on the wave vector c = 1

	type = 3;
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

const double* Laser::get_E_0() const {
	return e_0;
} 

const double* Laser::get_B_0() const {
	return b_0;
} 

void Laser::set_B_0(double new_B[3]){
	std::copy(new_B, new_B + 3, b_0);
}

// Calculate the electric and magnetic field considering the envelope of the laser
double* Laser::get_fields_envelope(double pos[3], double t){

	if (type != 3){throw std::runtime_error("Envelope not initalized");}

	double* fields = new double[6]; // E[3] x B[3] an array with the eletric and magnetic fields

	double xdk = pos[0] * k[0] + pos[1] * k[1] + pos[2] * k[2]; // k.x for future calculations

	double field_osc = cos(xdk - freq * t);// The field part responsible for oscilation

 	// Calculate the phase normalized to the length so that the function is 0 at phase = 0 and phase = L
	double phase = (xdk - env_freq * t)/length;

	double envelope = 0;

	if (phase < 0|| phase > M_PI){
		// Change the sin * sin to other functions or to a different function that has to be defined by a character probably
		envelope = sin(phase) * sin(phase); // Calculate the value of the phase | Envelope function
	} else {
		envelope = 0;
	}

	for (int i = 0; i < 3; i++){
		fields[i] = e_0[i] * field_osc * envelope;
	}
	for (int i = 3; i < 6; i++){
		fields[i] = b_0[i-3] * field_osc * envelope;
	}

	return fields;
}

double* Laser::get_fields(double* t, double pos[3], int* _type){

	if (_type){
		type = *_type;
	}

	if (type == 0){

		double* fields = new double[6];
		for (int i = 0; i < 3; ++i) fields[i] = e_0[i];
    	for (int i = 3; i < 6; ++i) fields[i] = b_0[i - 3];
    	return fields;

	} else if (type == 1){

		if (!t){throw std::runtime_error("No time given, but type = 1");}
		if (freq == 0){throw std::runtime_error("No frequency given, but type = 1");}

		double* fields = new double[6];

		double aux = cos(freq * (*t));

		for (int i = 0; i < 3; ++i) fields[i] = e_0[i] * aux;
    	for (int i = 3; i < 6; ++i) fields[i] = b_0[i - 3] * aux;
    	return fields; 

	} else if (type == 2){

		if (!t){throw std::runtime_error("No time given, but type = 2");}
		if (!pos){throw std::runtime_error("No position given, but type = 2");}
		if (freq == 0){throw std::runtime_error("No frequency given, but type = 2");}

		double* fields = new double[6]; // E[3] x B[3] an array with the eletric and magnetic fields

		double xdk = pos[0] * k[0] + pos[1] * k[1] + pos[2] * k[2]; // k.x for future calculations

		double field_osc = cos(xdk - freq * (*t));// The field part responsible for oscilation

		for (int i = 0; i < 3; i++) fields[i] = e_0[i] * field_osc;
		for (int i = 3; i < 6; i++) fields[i] = b_0[i - 3] * field_osc;
		

		return fields;

	} else if (type == 3){

		if (!t){throw std::runtime_error("No time given, but type = 3");}
		if (!pos){throw std::runtime_error("No position given, but type = 3");}
		if (freq == 0){throw std::runtime_error("No frequency given, but type = 3");}

		double* fields = new double[6]; // E[3] x B[3] an array with the eletric and magnetic fields

		double xdk = pos[0] * k[0] + pos[1] * k[1] + pos[2] * k[2]; // k.x for future calculations

		double field_osc = cos(xdk - freq * (*t));// The field part responsible for oscilation

 		// Calculate the phase normalized to the length so that the function is 0 at phase = 0 and phase = L
		double phase = (xdk - env_freq * (*t))/length;

		double envelope = 0;

		if (phase < 0|| phase > M_PI){
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