#include "functions.hpp"

// CONSTANTS (almost)

double MASS = 1; 
double CHARGE = 1;
double G = 1; // hum whats this one?
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

void boris(Particle*& particles, Laser* lasers, double time, double time_step, int n_of_particles,  int n_of_lasers){

	for (int i = 0; i < n_of_particles; ++i) {

		double E_field[3] = {0,0,0};
		double B_field[3] = {0,0,0};

		Particle& p = particles[i];

		const double* p_position = p.getPosition();
		const double* p_momentum = p.getMomentum();
		//const double* p_spin = p.getSpin();

		for (int j = 0; j < n_of_lasers; ++j){

			// Cycle to add and get the s Electric and Magnetic Fields

   			const double* laser_E = lasers[j].get_E_0();
   			const double* laser_B = lasers[j].get_B_0();

   			double* fields = lasers[j].get_fields(&time, p_position); // These might not be used

			for (int k = 0; k < 3; ++k) {
				E_field[k] = E_field[k] + fields[k];
				B_field[k] = B_field[k] + fields[k+3];
			}
		}

		// 1st E half-step
		double u_minus[3];
		for (int j = 0; j < 3; ++j) {
			u_minus[j] = p_momentum[j] + (CHARGE * time_step / (2.0 * MASS)) * E_field[j];
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

    }
}

// Creates the particles to be used in the simulation
void createParticles(Particle* particles, int particle_number){
	double pos[3] = {1,0,0};
	double mom[3] = {1,0,0};
    double spin[3] = {-1,0,1};

	for (int i = 0; i < particle_number; ++i) {
	    particles[i] = Particle(pos,mom,spin,i);  // Dynamically allocate objects
	}
}

// Creates a histogram based on the Diagnostics specified before
void PerformDiagnostics(std::vector<int>*& hist, Particle particle,  \
	std::string* params, double* bsize, double* bmax, double* bmin, int n_of_pars){

	if (n_of_pars == 0){
        throw std::runtime_error("Not Enough Parameters to perform Diagnostics!");
    } 

	// Perform the cycle to save the data to an unordered array with the counting bins
	
	for (int k = 0; k < n_of_pars; k++){
		double value = 0;
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

    	int index = std::round((value - bmin[k]) / bsize[k]);  // Compute bin index
    	hist[k][index]++;  // Increment corresponding bin
	}
}

// Setups the variables for the simulation and diagnostics
void setupInputVariable(std::ifstream& input_file, int& particle_n, double& timestep, double& totaltime, std::string*& params,
	 double*& binsize, double*& binmax, double*& binmin, double*& bin_n, int& n_par, Laser*& lasers, int& laser_number){

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

    int l_index = 1; // Laser index / amount - 1 

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

    	// Laser creation stuff
    	if ("L"+std::to_string(l_index) == it->first && !it->second.empty()){ 

    		int type = std::stoi(it->second);

    		// Check for initial eletric field 
    		if (values["E" + std::to_string(l_index)].empty()) {throw std::runtime_error("No Electric field initalized! INITALIZE IT :) ");}

    		if (type == 0){

    			if (!values["K" + std::to_string(l_index)].empty() && values["B" + std::to_string(l_index)].empty()){
    				std::cout << l_index<<"\n";
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

    				lasers[l_index - 1] = Laser(temp_E, temp_k, -l_index, 'k', 'y');

    			} else if (values["K" + std::to_string(l_index)].empty() && !values["B" + std::to_string(l_index)].empty()){

    				if (values["FREQ" + std::to_string(l_index)].empty()){
    					throw std::runtime_error("No frequency initalized for a B field time dependent field (type 1)");}

    				double temp_E[3] = {0}, temp_B[3] = {0};
    				double freq = std::stod(values["FREQ" + std::to_string(l_index)]);
    				// Parse the input file text to vectors for Eletric and Wave vector
    				parseVector(values["E" + std::to_string(l_index)], temp_E);
    				parseVector(values["B" + std::to_string(l_index)], temp_B);

    				lasers[l_index - 1] = Laser(temp_E, temp_B, -l_index, 'b', 'y', &freq);

    			} else {throw std::runtime_error("Initalized both K and B for a constant field (type 1) (or both are empty)");}

    		}else if (type == 2){

    			if (values["K" + std::to_string(l_index)].empty()){throw std::runtime_error("No wave vector (k), defined for a continuous laser (type 2)");}

    			double temp_E[3] = {0}, temp_k[3] = {0};

    			// Parse the input file text to vectors for Eletric and Wave vector
    			parseVector(values["E" + std::to_string(l_index)], temp_E);
    			parseVector(values["K" + std::to_string(l_index)], temp_k);

    			lasers[l_index - 1] = Laser(temp_E, temp_k, -l_index);

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

    			lasers[l_index - 1] = Laser(temp_E, temp_k, -l_index, env_l, env_freq);

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

    n_par -= 1;

    for (int i = 0; i < n_par; i++){

        bin_n[i] = (binmax[i] - binmin[i]) / binsize[i] + 1; // Calculate the size of the histogram

        if (bin_n[i] < 1){
        	throw std::runtime_error("Invalid Bin parameters for BIN_NUMBER= " +  std::to_string(i+1)); // Throw the error if the number of bins doesnt make sense
        }
    }
}

// Writes the wanted values of a given particle to a file
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

// Helper function to help parse vectors
void parseVector(const std::string& value, double vec[3]) {
    std::istringstream ss(value);
    char comma; // To skip the commas
    ss >> vec[0] >> comma >> vec[1] >> comma >> vec[2];
}