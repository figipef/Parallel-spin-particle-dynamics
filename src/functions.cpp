#include "functions.hpp"

int test(int a){
	return a + 123;
}

double* cross(double a[3], double b[3]) {

}

double gamma(double u[3]){

}

void boris(Particle*& particles, double time_step, double E_field[3], double B_field[3], int n_of_particles){

}

void createParticles(Particle* particles, int particle_number){
	double mom[3] = {1,1.5,0};
    double pos[3] = {1,1.5,0};
    double spin[3] = {0,1,0};

	for (int i = 0; i < particle_number; ++i) {
	    particles[i] = Particle(mom,pos,spin,i);  // Dynamically allocate objects
	}

}

void PerformDiagnostics(std::ifstream& output_file, Particle* particles,  \
	int particle_number, std::string par1, std::string par2, std::string par3, double b1_size, double b2_size, double b3_size, \
	double b1_max, double b2_max, double b3_max, double b1_min, double b2_min, double b3_min){

}