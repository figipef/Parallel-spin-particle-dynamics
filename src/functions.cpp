#include "functions.hpp"

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

void boris(Particle*& particles, double time_step, double E_field[3], double B_field[3], int n_of_particles){

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
