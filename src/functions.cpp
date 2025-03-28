#include "functions.hpp"
#include <memory>


// CONSTANTS (almost)

double MASS = 1; 
double CHARGE = 1;
double omega0 = 1;
double C = 1;
double G = 1; // hum whats this one? idk
double a = 0.00116; //anomalous magetic moment

// Calculate the cross product of two vectors
double* cross(double a[3], double b[3]) {
	double* c = new double[3];
	c[0] = a[1]*b[2] - a[2]*b[1];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];
	return c;
}

// Calculate the inner product of two vectors
double inner(double a[3], double b[3]){
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// Calculate the lorentz factor of a particle with generalized velocity u
double gamma(double u[3]){
	return sqrt(1. + inner(u, u));
}

//Calculate the spin precession vector, according to the BMT eq, all quantities, including gamma, should be in n+1/2 time (momentum time)
double* Omega(double gamma, double u[3], double E[3], double B[3]){
	double* Om = new double[3];
	double* ucE = cross(u, E); //u cross E
	double udB = inner(u, B); // u dot B
	for (int i = 0; i < 3; i++){
		Om[i] = ((a + 1./gamma)*(ucE[i] - gamma*B[i]) + a*udB*u[i]/(gamma + 1.))/gamma;
	}
	return Om;
}

double* f_RR(double gamma, double u[3], double E[3], double B[3]){
	double dotB_u = inner(B, u);
	double dotE_u = inner(E,u);
 	double b2 = inner(B, B);
	double* crossB_E = cross(B,E);
	double sigma0 = (2*omega0*CHARGE*CHARGE)/(3*MASS*C*C*C);

	// defining F^2u
	double F2_u[3]; 
	for (int i = 0; i<3; i++){
		F2_u[i] = dotE_u * E[i] + crossB_E[i] * u[i] + B[i] * dotB_u - b2 * u[i];
	}  

	// (u | F^2 u)
	double dotu_F2u = -inner(u,F2_u);
	// (u | F^2 u)u
	double spatial_result[3];
	for (int i = 0; i<3; i++){
		spatial_result[i] = F2_u[i] - dotu_F2u * u[i];
	} 

	double factor = (sigma0 * CHARGE * CHARGE) / (gamma * MASS);
	double* result = new double[3];
    for (int i = 0; i < 3; i++) {
        result[i] = factor * spatial_result[i];
    }
    
    return result;

}

// Creates a histogram based on the Diagnostics specified before
void PerformDiagnostics(Histogram& hist, Particle particle,  \
	std::string* params, double* bsize, double* bmax, double* bmin, int n_of_pars){

	if (n_of_pars == 0){
        throw std::runtime_error("Not Enough Parameters to perform Diagnostics!");
    } 

	// Perform the cycle to save the data to an unordered array with the counting bins
	
	if (n_of_pars == 1) {

		double value = 0;
		int i = params[0][1] - '0'; // Get the particle index from the parameter character
		switch(params[0][0]){
			case 'p':
				value = particle.getPosition()[i - 1];
				break;
			case 'm':
				value = particle.getMomentum()[i - 1];
				break;
			case 's':
				value = particle.getSpin()[i - 1];
				break;
		} 

		if (!(value < bmin[0] || value > bmax[0])) { // Ignore out-of-range values

			int index = std::round((value - bmin[0]) / bsize[0]); // Compute bin index
			hist.vector1D[index]++; // Increment corresponding bin
		}  

	} else if (n_of_pars == 2) {

		double value1 = 0;
		int i = params[0][1] - '0'; // Get the particle index from the parameter character
		switch(params[0][0]){
			case 'p':
				value1 = particle.getPosition()[i - 1];
				break;
			case 'm':
				value1 = particle.getMomentum()[i - 1];
				break;
			case 's':
				value1 = particle.getSpin()[i - 1];
				break;
		} 

		double value2 = 0;
		i = params[1][1] - '0'; // Get the particle index from the parameter character
		switch(params[1][0]){
			case 'p':
				value2 = particle.getPosition()[i - 1];
				break;
			case 'm':
				value2 = particle.getMomentum()[i - 1];
				break;
			case 's':
				value2 = particle.getSpin()[i - 1];
				break;
		}

		int index1;
		int index2;
		bool i1 = false;
		bool i2 = false;

		if (!(value1 < bmin[0] || value1 > bmax[0])) {index1 = std::round((value1 - bmin[0]) / bsize[0]); i1 = true;}  // Ignore out-of-range values and Compute bin index
		if (!(value2 < bmin[1] || value2 > bmax[1])) {index2 = std::round((value2 - bmin[1]) / bsize[1]); i2 = true;}  // Ignore out-of-range values and Compute bin index

		if (i1 && i2){
			hist.matrix2D[index2][index1]++;  // Increment corresponding bin
		}

	} else {

		throw std::runtime_error("PerformDiagnostics: Invalid number of params");
	}
}

void boris(Particle*& particles, Laser* lasers, double time, double time_step, int n_of_particles, int n_of_lasers, int RR, Histogram* hist, DiagnosticParameters* diag_params){

	for (int i = 0; i < n_of_particles; ++i) {

		double E_field[3] = {0,0,0};
		double B_field[3] = {0,0,0};

		Particle& p = particles[i];

		const double* p_position = p.getPosition();
		const double* p_momentum = p.getMomentum();
		//const double* p_spin = p.getSpin();

		for (int j = 0; j < n_of_lasers; ++j){

			// Cycle to add and get the s Electric and Magnetic Fields

   			//const double* laser_E = lasers[j].get_E_0();
   			//const double* laser_B = lasers[j].get_B_0();

   			double* fields = lasers[j].get_fields(&time, p_position); // These might not be used
   			//std::cout << "Pos: "<<p_position[0]<<" "<<p_position[1]<<" "<<p_position[2]<<"\n";

			for (int k = 0; k < 3; ++k) {
				E_field[k] = E_field[k] + fields[k];
				B_field[k] = B_field[k] + fields[k+3];
				//std::cout << "E " << E_field[k] << " B "<< B_field[k] <<"\n";
			}
		}

		double p_momentum_copy[3];
		p_momentum_copy[0] = p_momentum[0];
		p_momentum_copy[1] = p_momentum[1];
		p_momentum_copy[2] = p_momentum[2];

		double* f_RR_Value = f_RR(gamma(p_momentum_copy),p_momentum_copy,E_field,B_field);
	
		// 1st E half-step
		double u_minus[3];
	
		for (int j = 0; j < 3; ++j) {
			u_minus[j] = p_momentum[j] + (CHARGE * time_step / (2.0 * MASS)) * E_field[j] + RR * (time_step / 2.0)*f_RR_Value[j];
		};

        // auxiliary vector t
        double t_vec[3];
        double factor = CHARGE * time_step / (2.0 * MASS * gamma(u_minus));
        for (int j = 0; j < 3; ++j) {
            t_vec[j] = factor * B_field[j];
        }
        double t2 = inner(t_vec,t_vec);
        
        // auxiliary vector s 
        double s_vec[3];
        for (int j = 0; j < 3; ++j) {
            s_vec[j] = 2.0 * t_vec[j] / (1.0 + t2);
        }
        
        // B rotation
        double u_prime[3];
		double* uct = cross(u_minus, t_vec);
		for (int j = 0; j < 3; ++j) {
			u_prime[j] = u_minus[j] + uct[j];
		}
        
		double u_plus[3];
		double* ucs = cross(u_prime, s_vec);
		for (int j = 0; j < 3; ++j) {
			u_plus[j] = u_minus[j] + ucs[j];
		}
        
        // 2nd E half-step
        double u_next[3];
        for (int j = 0; j < 3; ++j) {
			u_next[j] = u_plus[j] + (CHARGE * time_step / (2.0 * MASS)) * E_field[j] + RR * (time_step / 2.0)*f_RR_Value[j];
		}
        
        // Update the particle's momentum
		p.setMomentum(u_next);

		// Now we update the particle position with the new momentum
		double new_gamma = gamma(u_next);

		double new_v[3];
		for (int j= 0; j <3; j++){
			new_v[j] = u_next[j] / (new_gamma * MASS);  
		}

		double new_pos[3];
		for (int j= 0; j <3; j++){
			new_pos[j] = p_position[j] + time_step * new_v[j]; 
		}

		particles[i].setPosition(new_pos);

		//Spin update

		// reusing auxiliary vectors t and s
        factor = time_step / 2.0;
		double* Om = Omega(new_gamma, u_next, E_field, B_field);
        for (int j = 0; j < 3; ++j) {
            t_vec[j] = factor * Om[j];
        }
        t2 = inner(t_vec,t_vec);

		for (int j = 0; j < 3; ++j) {
            s_vec[j] = 2.0 * t_vec[j] / (1.0 + t2);
        }

		double* spn = new double[3];  // Dynamically allocate memory
		std::copy(p.getSpin(), p.getSpin() + 3, spn);
		double spn_prime[3];
		double* sct = cross(spn, t_vec);
		for (int j = 0; j < 3; ++j) {
			spn_prime[j] = spn[j] + sct[j];
		}
        
		double spn_next[3];
		double* scs = cross(spn_prime, s_vec);
		for (int j = 0; j < 3; ++j) {
			spn_next[j] = spn[j] + scs[j];
		}
		p.setSpin(spn_next);

		// If the histogram is passed, perform the Diagnostics
		if(hist != nullptr && diag_params != nullptr){

			PerformDiagnostics(*hist, p, diag_params -> params, \
				diag_params -> bsize, diag_params -> bmax, diag_params -> bmin, diag_params -> n_of_pars);
		}
    }
}

// Creates the particles to be used in the simulation
void createParticles(Particle* particles, int particle_number, std::string* types, double* dist_sizes, double* momentum_dir, double* spin_dir, Laser* lasers, int n_of_lasers){ 

	// types are defined by "uni" and "gau" || probably eventually change this to 0s && 1s to run faster
	// CHANGE THIS NAME PLS dist_sizes is just the length of the box for uniform and a std for gaussian

	//Calculate the initial position according to the lasers

	int largest_laser_index = -1;
	double largest_laser_length = 0;
	double pos_shift[3] = {0,0,0};

	for (int i = 0; i < n_of_lasers; i++){

		if (lasers[i].type == 3 && lasers[i].length > largest_laser_length) {

			largest_laser_length = lasers[i].length;
			largest_laser_index = i;
		}
	}

	if (largest_laser_index != 1){ // Adjust the initial positon according to the size of the laser
		Laser l = lasers[largest_laser_index];

		double abs_k = std::sqrt(l.k[0] * l.k[0] + l.k[1] * l.k[1] + l.k[2] + l.k[2]);
		double mult = 0; // top differentiate from uniform and gaussian 

		if (types[0] == "0"){
			mult = dist_sizes[0];
		} else {
			mult = dist_sizes[0] * 3;
		}

		// To copy the sign for k in case of negative k's
		pos_shift[0] = l.k[0]/abs_k * l.length + (l.k[0] == 0 ? 0.0 : (std::copysign(1.0, l.k[0]) * mult));
		pos_shift[1] = l.k[1]/abs_k * l.length + (l.k[1] == 0 ? 0.0 : (std::copysign(1.0, l.k[1]) * mult)); 
		pos_shift[2] = l.k[2]/abs_k * l.length + (l.k[2] == 0 ? 0.0 : (std::copysign(1.0, l.k[2]) * mult));
		
	}

    double spin[3] = {0,0,0};

    std::random_device rd;   // Random seed generator
    std::mt19937 gen(rd());  // Mersenne Twister PRNG seeded with rd()

    // Declare distributions for future assignment
    std::unique_ptr<std::uniform_real_distribution<double>> space_uniform;
    std::unique_ptr<std::normal_distribution<double>> space_gaussian;
    std::unique_ptr<std::uniform_real_distribution<double>> momentum_uniform;
    std::unique_ptr<std::normal_distribution<double>> momentum_gaussian;
	std::unique_ptr<std::uniform_real_distribution<double>> unif_1to0_dist; //standard 0 to 1 uniform dist im using for spin

    if (types[0] == "0"){

    	space_uniform = std::make_unique<std::uniform_real_distribution<double>>(-dist_sizes[0], dist_sizes[0]);
    
    } else {

    	space_gaussian = std::make_unique<std::normal_distribution<double>>(0.0, dist_sizes[0]);
    }

    if (types[1] == "0") {

    	momentum_uniform = std::make_unique<std::uniform_real_distribution<double>>(-dist_sizes[1], dist_sizes[1]);

    } else {

    	momentum_gaussian = std::make_unique<std::normal_distribution<double>>(0.0, dist_sizes[1]);

    }

	for (int i = 0; i < particle_number; ++i) {

		double pos[3] = {0,0,0};
		double mom[3] = {0,0,0};


		// Assign the position values
		if (types[0] == "0"){

			pos[0] = (*space_uniform)(gen);
			pos[1] = (*space_uniform)(gen);
			pos[2] = (*space_uniform)(gen);

		} else if (types[0] == "1") { // For the Gaussian case, don't surpass 3 std's

			do { pos[0] = (*space_gaussian)(gen); } while (std::abs(pos[0]) > 3 * dist_sizes[0]);
			do { pos[1] = (*space_gaussian)(gen); } while (std::abs(pos[1]) > 3 * dist_sizes[0]);
			do { pos[2] = (*space_gaussian)(gen); } while (std::abs(pos[2]) > 3 * dist_sizes[0]);

		} else if (types[0] == "2") {
			throw std::runtime_error("Still undefined");
		} else { throw std::runtime_error("Invalid position type in particle creator"); }

		// Assign the momentum values
		if (types[1] == "0"){

			mom[0] = (*momentum_uniform)(gen);
			mom[1] = (*momentum_uniform)(gen);
			mom[2] = (*momentum_uniform)(gen);

		} else if (types[1] == "1"){ // For the Gaussian case, don't surpass 3 std's

			do { mom[0] = (*momentum_gaussian)(gen); } while (std::abs(mom[0]) > 3 * dist_sizes[1]);
			do { mom[1] = (*momentum_gaussian)(gen); } while (std::abs(mom[1]) > 3 * dist_sizes[1]);
			do { mom[2] = (*momentum_gaussian)(gen); } while (std::abs(mom[2]) > 3 * dist_sizes[1]);

		} else if (types[1] == "2"){

			mom[0] = momentum_dir[0];
			mom[1] = momentum_dir[1];
			mom[2] = momentum_dir[2];

		} else { throw std::runtime_error("Invalid position type in particle creator"); }

		pos[0] = pos[0] + pos_shift[0];
		pos[1] = pos[1] + pos_shift[1];
		pos[2] = pos[2] + pos_shift[2];

		if (types[2] == "0"){

			unif_1to0_dist = std::make_unique<std::uniform_real_distribution<double>>(0., 1.);

			double phi = 2.*M_PI*(*unif_1to0_dist)(gen);
			double cos_th = 2.*(*unif_1to0_dist)(gen) - 1.;
			double sin_th = sqrt(1.- cos_th*cos_th);

			spin[0] = sin_th*cos(phi);
			spin[1] = sin_th*sin(phi);
			spin[2] = cos_th;
		
		} else if (types[2] == "1"){

			double Norm = sqrt(inner(spin_dir, spin_dir)); //sqrt(spin_dir[0]*spin_dir[0] + spin_dir[1]*spin_dir[1] + spin_dir[2]*spin_dir[2]);
			spin[0] = spin_dir[0] / Norm;
			spin[1] = spin_dir[1] / Norm;
			spin[2] = spin_dir[2] / Norm;
        
		} else if (types[2] == "2"){
			//double* s_try = new double[3];
			bool accept = false;
			unif_1to0_dist = std::make_unique<std::uniform_real_distribution<double>>(0., 1.);
			while (!accept){
				
				double phi = 2.*M_PI*(*unif_1to0_dist)(gen);
				double cos_th = 2.*(*unif_1to0_dist)(gen) - 1.;
				double sin_th = sqrt(1.- cos_th*cos_th);
	
				spin[0] = sin_th*cos(phi);
				spin[1] = sin_th*sin(phi);
				spin[2] = cos_th;

				double Norm = sqrt(inner(spin_dir, spin_dir));
				for (int d=0; d<3; d++){
					spin_dir[d] = spin_dir[d]/Norm;
				}
				

				double prob = (*unif_1to0_dist)(gen);
				double val = exp(dist_sizes[2]*(inner(spin_dir, spin)))*dist_sizes[2]/(4*M_PI*sinh(dist_sizes[2]));
				if (prob <= val){accept = true;}
			}

		}

	    particles[i] = Particle(pos,mom,spin,i);  // Create the particle to the array

	}
}

void FieldDiagWritter(double& dt, int& iter, double*& fieldiag, Laser*& lasers, int& laser_number){
	double dx, min, max, x;
	double t;
	int dir;
	double pos[3] = {0.,0.,0.};

	if (fieldiag[0] == 1.){
		dx = fieldiag[2]; 
		min = fieldiag[3];
		max = fieldiag[4];
		dir = fieldiag[5];
		std::ofstream e_field_1("../output/e_field1.txt");
		std::ofstream e_field_2("../output/e_field2.txt");
		std::ofstream e_field_3("../output/e_field3.txt");
		std::cout<<"look here idiot\n";
		std::cout << "E = [ " <<lasers[0].get_E_0()[0]<<", "<<lasers[0].get_E_0()[1]<<", "<<lasers[0].get_E_0()[2]<< " ]\n";
        std::cout << "B = [ " <<lasers[0].get_B_0()[0]<<", "<<lasers[0].get_B_0()[1]<<", "<<lasers[0].get_B_0()[2]<< " ]\n";
		double t1 = 0.;
		double pos1[3] = {2.,0,0};
		for (t = 0.; t < iter*dt; t = t + dt){
			for (x = min; x <= max; x = x + dx){
				pos[dir-1] = x;
				e_field_1 << lasers[0].get_fields(&t,pos)[0] <<" ";
				e_field_2 << lasers[0].get_fields(&t,pos)[1] <<" ";
				e_field_3 << lasers[0].get_fields(&t,pos)[2] <<" ";
			}
			e_field_1 << "\n";
			e_field_2 << "\n";
			e_field_3 << "\n";
		}
	}

	if (fieldiag[1] == 1.){
		dx = fieldiag[2]; 
		min = fieldiag[3];
		max = fieldiag[4];
		dir = fieldiag[5];
		std::ofstream b_field_1("../output/b_field1.txt");
		std::ofstream b_field_2("../output/b_field2.txt");
		std::ofstream b_field_3("../output/b_field3.txt");

		for (t = 0.; t < iter*dt; t = t + dt){
			for (x = min; x <= max; x = x + dx){
				pos[dir-1] = x;
				b_field_1 << lasers[0].get_fields(&t,pos)[3] <<" ";
				b_field_2 << lasers[0].get_fields(&t,pos)[4] <<" ";
				b_field_3 << lasers[0].get_fields(&t,pos)[5] <<" ";
			}
			b_field_1 << "\n";
			b_field_2 << "\n";
			b_field_3 << "\n";
		}
	}
}

// Setups the variables for the simulation and diagnostics
void setupInputVariable(std::ifstream& input_file, int& particle_n, std::string*& dist_types, double*& dist_sizes, double*& momentum_dir, double*& spin_dir, double& timestep, double& totaltime, int& step_diag, std::string*& params,
	 double*& binsize, double*& binmax, double*& binmin, int*& bin_n, int& n_par, double*& fieldiag, Laser*& lasers, int& laser_number, bool& RR){

	/*
	Function that takes both simulation and diagnostics input parameters and updates them through the input_file
	*/

    std::unordered_map<std::string, std::string> values;
    std::string line;
    
    // Extract the values line by line from the input file and store
    while (std::getline(input_file, line)) {
        std::istringstream iss(line);
        std::string key, value, unit;
        
        if (std::getline(iss, key, '=')) {
            iss >> value >> unit; // Extract value and unit
            values[key] = value;  // Store value in the map
        }
    }
    input_file.close();

    // Save the input values to usable variables

	RR = std::stoi(values["RADIATION_REACTION"]);
    particle_n = std::stoi(values["NUMBER_OF_PARTICLES"]);
    timestep = std::stod(values["TIME_STEP"]);
    totaltime = std::stod(values["TOTAL_TIME"]);

    step_diag = std::stoi(values["STEPS_DIAG"]);

    // Assign particle creation distribution parameters
    dist_types[0] = values["POSITION_DIST_TYPE"];
    dist_types[1] = values["MOMENTUM_DIST_TYPE"];
	dist_types[2] = values["SPIN_DIST_TYPE"];

	if (dist_types[1] == "2"){
		if(values["MOMENTUM_DIST_TYPE"].empty()){{throw std::runtime_error("Preferred Momentum Direction Distribution Selected But No Direction Indicated! Please Initialize It Correctly in Input");}}
		parseVector(values["MOMENTUM_PREF_DIR"], momentum_dir);
	}

	if (dist_types[2] == "1" || dist_types[2] == "2"){
		if(values["SPIN_PREF_DIR"].empty()){{throw std::runtime_error("Preferred Spin Direction Distribution Selected But No Direction Indicated! Please Initialize It Correctly in Input");}}
		parseVector(values["SPIN_PREF_DIR"], spin_dir);
	}

    dist_sizes[0] = std::stod(values["POSITION_DIST_SIZE"]);
    dist_sizes[1] = std::stod(values["MOMENTUM_DIST_SIZE"]);

	if (dist_types[2] == "2"){
		if(values["SPIN_DIST_SIZE"].empty()){{throw std::runtime_error("Von Mises-Fisher Distribution Selected But No Size Parameter Given! Please Initialize It Correctly in Input");}}
		dist_sizes[2] = std::stod(values["SPIN_DIST_SIZE"]);
	}

    // Initalize counters for the bin parameters

    n_par = 1;
    int nbs = 1;
    int bM = 1;
    int bm = 1;

    int l_index = 1; // Laser index / amount - 1 

    bool show_b = true;
    bool show_e = true;
    bool field_bin = true;
    bool field_bin_min = true;
    bool field_bin_max = true;
    bool field_bin_dir = true;

    auto it = values.begin(); // Iterator for the map

    while (it != values.end()) {

    	if ("PAR"+std::to_string(n_par) == it->first && !it->second.empty()){

    		// Checkers to accept correct param definition
    		if (it->second.length() != 2){
    			throw std::runtime_error("Invalid Parameter for Parameter= " +  std::to_string(n_par) + "; parameter not the right size");
    		}

    		if (it->second[0] != 's' && it->second[0] != 'm' && it->second[0] != 'p'){
    			throw std::runtime_error("Invalid Parameter for Parameter= " +  std::to_string(n_par) + "; not a valid character at position 0");
    		}

    		if (!std::isdigit(it->second[1])){
    			throw std::runtime_error("Invalid Parameter for Parameter= " +  std::to_string(n_par) + "; not a digit at position 1");
    		}

    		params[n_par - 1] = it->second; // Save it to the array
    		n_par += 1; // Counter update
    		it = values.begin(); // Restart Iterator
    	}

    	if ("BIN_NUMBER_"+std::to_string(nbs) == it->first && !it->second.empty()){

    		bin_n[nbs - 1] =std::stoi(it->second);// Save it to the array

    		if (bin_n[nbs - 1] < 1){ // Check to see if the number of bins is valid
    			throw std::runtime_error("Invalid Number of Bins <1 ");
    		}

    		nbs += 1;// Counter update
    		it = values.begin();// Restart Iterator
    	}

    	if ("BIN_MAX_"+std::to_string(bM) == it->first && !it->second.empty()){
    		
    		binmax[bM - 1] =std::stod(it->second); // Save it to the array
    		bM += 1; // Counter update
    		it = values.begin(); // Restart Iterator
    	}

    	if ("BIN_MIN_"+std::to_string(bm) == it->first && !it->second.empty()){

    		binmin[bm - 1] =std::stod(it->second); // Save it to the array
    		bm += 1; // Counter update
    		it = values.begin(); // Restart Iterator
    	}
		
		if ("SHOW_E" == it->first && !it->second.empty() && show_e){

    		if ("1" == it->second ){fieldiag[0] = 1.;}
			else {fieldiag[0] = 0.;}
			show_e = false;
			it = values.begin(); // Restart Iterator

    	}

		if ("SHOW_B" == it->first && !it->second.empty() && show_b){

    		if ("1" == it->second ){fieldiag[1] = 1.;}
			else {fieldiag[1] = 0.;}
			show_b = false;
			it = values.begin(); // Restart Iterator
			
    	}

		if ("FIELD_BIN" == it->first && field_bin && ( fieldiag[0] == 1. || fieldiag[1] == 1)){

    		fieldiag[2] = std::stod(it->second);
    		field_bin = false;
    		it = values.begin(); // Restart Iterator

    	}

		if ("FIELD_BIN_MIN" == it->first && field_bin_min && ( fieldiag[0] == 1. || fieldiag[1] == 1)){

    		fieldiag[3] = std::stod(it->second);
    		field_bin_min = false;
    		it = values.begin(); // Restart Iterator
    		
    	}

		if ("FIELD_BIN_MAX" == it->first && field_bin_max && ( fieldiag[0] == 1. || fieldiag[1] == 1) ){

    		fieldiag[4] = std::stod(it->second);
    		field_bin_max = false;
    		it = values.begin(); // Restart Iterator
    		
    	}

		if ("FIELD_BIN_DIR" == it->first && field_bin_dir && ( fieldiag[0] == 1. || fieldiag[1] == 1)){

    		fieldiag[5] = std::stoi(it->second);
    		field_bin_dir = false;
    		it = values.begin(); // Restart Iterator
    		
    	}	

    	// Laser creation stuff
    	if ("L"+std::to_string(l_index) == it->first && !it->second.empty()){ 

    		int type = std::stoi(it->second);
    		int ext_phase = 0;

    		// Set the external phase of the laser
    		if (!values["PHASE" + std::to_string(l_index)].empty()){
    			ext_phase = std::stoi(values["PHASE" + std::to_string(l_index)]);

    		}

    		// Check for initial eletric field 
    		if (values["E" + std::to_string(l_index)].empty()) {throw std::runtime_error("No Electric field initalized! INITALIZE IT :) ");}

    		if (type == 0){

    			if (!values["K" + std::to_string(l_index)].empty() && values["B" + std::to_string(l_index)].empty()){
    				double temp_E[3] = {0}, temp_k[3] = {0};

    				// Parse the input file text to vectors for Eletric and Wave vector
    				parseVector(values["E" + std::to_string(l_index)], temp_E);
    				parseVector(values["K" + std::to_string(l_index)], temp_k);

    				lasers[l_index - 1] = Laser(temp_E, temp_k, -l_index, 'k', 'n');

    			} else if (values["K" + std::to_string(l_index)].empty() && !values["B" + std::to_string(l_index)].empty()){

    				double temp_E[3] = {0}, temp_B[3] = {0};

    				// Parse the input file text to vectors for Eletric and Wave vector
    				parseVector(values["E" + std::to_string(l_index)], temp_E);
    				parseVector(values["B" + std::to_string(l_index)], temp_B);

    				lasers[l_index - 1] = Laser(temp_E, temp_B, -l_index, 'b', 'n');

    			} else {throw std::runtime_error("Initalized both K and B for a constant field (type 0)");}

    		}else if (type == 1){

    			if (!values["K" + std::to_string(l_index)].empty() && values["B" + std::to_string(l_index)].empty()){

    				double temp_E[3] = {0}, temp_k[3] = {0};

    				// Parse the input file text to vectors for Eletric and Wave vector
    				parseVector(values["E" + std::to_string(l_index)], temp_E);
    				parseVector(values["K" + std::to_string(l_index)], temp_k);

    				lasers[l_index - 1] = Laser(temp_E, temp_k, -l_index, 'k', 'y', ext_phase);

    			} else if (values["K" + std::to_string(l_index)].empty() && !values["B" + std::to_string(l_index)].empty()){

    				if (values["FREQ" + std::to_string(l_index)].empty()){
    					throw std::runtime_error("No frequency initalized for a B field time dependent field (type 1)");}

    				double temp_E[3] = {0}, temp_B[3] = {0};
    				double freq = std::stod(values["FREQ" + std::to_string(l_index)]);
    				// Parse the input file text to vectors for Eletric and Wave vector
    				parseVector(values["E" + std::to_string(l_index)], temp_E);
    				parseVector(values["B" + std::to_string(l_index)], temp_B);

    				lasers[l_index - 1] = Laser(temp_E, temp_B, -l_index, 'b', 'y', ext_phase, &freq);

    			} else {throw std::runtime_error("Initalized both K and B for a constant field (type 1) (or both are empty)");}

    		}else if (type == 2){

    			if (values["K" + std::to_string(l_index)].empty()){throw std::runtime_error("No wave vector (k), defined for a continuous laser (type 2)");}

    			double temp_E[3] = {0}, temp_k[3] = {0};

    			// Parse the input file text to vectors for Eletric and Wave vector
    			parseVector(values["E" + std::to_string(l_index)], temp_E);
    			parseVector(values["K" + std::to_string(l_index)], temp_k);

    			lasers[l_index - 1] = Laser(temp_E, temp_k, -l_index, ext_phase);

    		}else if (type == 3){

    			if (values["K" + std::to_string(l_index)].empty()){
    				throw std::runtime_error("No wave vector (k) defined for a packet laser (type 3)");}

    			if (values["ENV_FREQ" + std::to_string(l_index)].empty()){
    				throw std::runtime_error("No envelope frequency defined for a packet laser (type 3)");}
    				
    			if (values["ENV_L" + std::to_string(l_index)].empty()){
    				throw std::runtime_error("No envelope length defined for a packet laser (type 3)");}	

    			double env_l = std::stod(values["ENV_L" + std::to_string(l_index)]);
    			double env_freq = std::stod(values["ENV_FREQ" + std::to_string(l_index)]);

    			double temp_E[3] = {0}, temp_k[3] = {0};

    			// Parse the input file text to vectors for Eletric and Wave vector
    			parseVector(values["E" + std::to_string(l_index)], temp_E);
    			parseVector(values["K" + std::to_string(l_index)], temp_k);

    			lasers[l_index - 1] = Laser(temp_E, temp_k, -l_index, env_l, env_freq, ext_phase);

    		} else {throw std::runtime_error("Invalid Laser type");}

    		l_index++;
    	}

    	it ++;
    }

    laser_number = l_index - 1; 

    // Check if the number of variables is the same for every bin parameter and parameters
    if (n_par != bm || n_par != bM || n_par != nbs){
    	throw std::runtime_error("Number of Input parameters does not match the bin parameters");
    }

    if (n_par > 3){
    	throw std::runtime_error("Number of diagnostic parameters larger than 2");
    }

    n_par -= 1;

    for (int i = 0; i < n_par; i++){

    	binsize[i] = (binmax[i] - binmin[i]) / (bin_n[i] - 1 + 1e-308); // calculate the size of each bin

        //bin_n[i] = std::round((binmax[i] - binmin[i]) / binsize[i] + 1); // Calculate the size of the histogram
    	//std::cout << binsize[i];
        if (binsize[i] < 0){
        	throw std::runtime_error("Invalid Bin parameters for BIN_NUMBER= " +  std::to_string(i+1)); // Throw the error if the number of bins doesnt make sense
        }
    }
}
