#include "functions.hpp"

// CONSTANTS (almost)

double MASS = 1; 
double CHARGE = 1;
double G = 1;
double a = 0.00116; //anomalous magetic moment

int test(int a){
	return a + 123;
}

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

void boris(Particle*& particles, Laser* lasers, double time_step, int n_of_particles,  int n_of_lasers){
 	// NECESSÁRIO ADICIONAR TEMPO!!!! PARA DEPOIS CALCULAR OS CAMPOS DEPENDENTES DO TEMPO 

	for (int i = 0; i < n_of_particles; ++i) {

		double E_field[3] = {0,0,0};
		double B_field[3] = {0,0,0};

		for (int j = 0; j < n_of_lasers; ++j){

			// Cycle to add and get the s Electric and Magnetic Fields

   			const double* laser_E = lasers[j].get_E_0();
   			const double* laser_B = lasers[j].get_B_0();

			for (int k = 0; k < 3; ++k) {
				E_field[k] = E_field[k] + laser_E[k];
				B_field[k] = B_field[k] + laser_B[k];
			}
		}

		Particle& p = particles[i];

		// 1st E half-step
		double u_minus[3];
		for (int j = 0; j < 3; ++j) {
			u_minus[j] = p.getMomentum()[j] + (CHARGE * time_step / (2.0 * MASS)) * E_field[j];
		}

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
            u_next[j] = u_plus[j] + (CHARGE * time_step / (2.0 * MASS)) * E_field[j];
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
			new_pos[j] = particles[i].getPosition()[j] + time_step * new_v[j]; 
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

    }
}

void createParticles(Particle* particles, int particle_number){
	double pos[3] = {1,1.5,0};
	double mom[3] = {1,1.5,0};
    double spin[3] = {-1,0,1};

	for (int i = 0; i < particle_number; ++i) {
	    particles[i] = Particle(mom,pos,spin,i);  // Dynamically allocate objects
	}
}

void PerformDiagnostics(std::vector<int>*& hist, Particle* particles,  \
	int particle_number, std::string par1, std::string par2, std::string par3, double b1_size, double b2_size, double b3_size, \
	double b1_max, double b2_max, double b3_max, double b1_min, double b2_min, double b3_min){

	// Initialize the arrays for relevant quantities
	char *p = new char[3];
	int *i = new int[3];
	double *bmax = new double[3];
	double *bmin = new double[3];
	double *bsize = new double[3];

	int n_of_pars = 0;

	// Check if the parameters have values and stores the number of parameters needed
	if (!par1.empty()){
		n_of_pars += 1;

		p[0] = par1[0];
		i[0] = par1[1];
		bmax[0] = b1_max;
		bmin[0] = b1_min;
		bsize[0] = b1_size;
	}

	if (!par2.empty()){
		n_of_pars += 1;

		p[1] = par2[0];
		i[1] = par2[1];
		bmax[1] = b2_max;
		bmin[1] = b2_min;
		bsize[1] = b2_size;

	}

	if (!par3.empty()){
		n_of_pars += 1;

		p[2] = par3[0];
		i[2] = par3[1];
		bmax[2] = b3_max;
		bmin[2] = b3_min;
		bsize[2] = b3_size;
	}

	// Perform the cycle to save the data to an unordered array with the counting bins
	for (int j = 0; j < particle_number; j++){
		for (int k = 0; k < n_of_pars; k++){
			double value;
			switch(p[k]){
				case 'p':
					value = particles[j].getPosition()[i[k]];
					break;
				case 'm':
					value = particles[j].getMomentum()[i[k]];
					break;
				case 's':
					value = particles[j].getSpin()[i[k]];
					break;
			}

			if (value < bmin[k] || value > bmax[k]) break;  // Ignore out-of-range values

    		int index = (value - bmin[k]) / bsize[k];  // Compute bin index
    		hist[k][index]++;  // Increment corresponding bin
		}
	}
}

void PerformDiagnostics(std::vector<int>*& hist, Particle particle,  \
	std::string* params, double* bsize, double* bmax, double* bmin, int n_of_pars){

	if (n_of_pars == 0){
		
        throw std::runtime_error("Not Enough Parameters to perform Diagnostics!");
    } 

	// Perform the cycle to save the data to an unordered array with the counting bins
	
	for (int k = 0; k < n_of_pars; k++){
		double value;
		int i = params[k][1] - '0'; // Get the particle index from the parameter character
		switch(params[k][0]){
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

		if (value < bmin[k] || value > bmax[k]) break;  // Ignore out-of-range values

    	int index = (value - bmin[k]) / bsize[k];  // Compute bin index
    	hist[k][index]++;  // Increment corresponding bin
	}
}

void setupInputVariable(std::ifstream& input_file, int& particle_n, double& timestep, double& totaltime, std::string*& params, double*& binsize, double*& binmax, double*& binmin, double*& bin_n, int& n_par){

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

    particle_n = std::stoi(values["NUMBER_OF_PARTICLES"]);
    timestep = std::stod(values["TIME_STEP"]);
    totaltime = std::stod(values["TOTAL_TIME"]);

    // Initalize counters for the bin parameters

    n_par = 1;
    int nbs = 1;
    int bM = 1;
    int bm = 1;

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

    	if ("BIN_SIZE_"+std::to_string(nbs) == it->first && !it->second.empty()){

    		binsize[nbs - 1] =std::stod(it->second);// Save it to the array
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
    	it ++;
    }

    // Check if the number of variables is the same for every bin parameter and parameters
    if (n_par != bm || n_par != bM || n_par != nbs){
    	throw std::runtime_error("Number of Input parameters does not match the bin parameters");
    }

    n_par -= 1;

    for (int i = 0; i < n_par; i++){

        bin_n[i] = (binmax[i] - binmin[i]) / binsize[i] + 1; // Calculate the size of the histogram

        if (bin_n[i] < 1){
        	throw std::runtime_error("Invalid Bin parameters for BIN_NUMBER= " +  std::to_string(i+1)); // Throw the error if the number of bins doesnt make sense
        }
    }
}

void writeToFile(std::ofstream& file, const Particle& p, char a){
	const double* data = nullptr;

	switch(a){

		case 'p':

			data = p.getPosition();
			break;
				
		case 'm':
			data = p.getMomentum();
				
			break;

		case 's':

			data = p.getSpin();
			break;

		default:
	        file << "Invalid option" << std::endl;
            break;
		}

	// Check for nullptr before accessing data
	if (data) {

	    for (int i = 0; i < 3; i++) {
	        file << data[i] << " ";
	    }
	    file << std::endl;

	} else {
	    file << "Error: Null data pointer" << std::endl;
	}	
}