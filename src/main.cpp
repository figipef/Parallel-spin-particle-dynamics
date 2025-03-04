#include "particle.hpp"
#include "functions.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <vector>

// CONSTANTS (almost)
double MASS = 1; 
double CHARGE = 1;
double G = 1;

int main() {

    // ============================
    //   VARIABLE INTIALIZATION
    // ============================

	std::ifstream input_file("input.txt"); // Open the input file

    if (!input_file) {
        throw std::runtime_error("Could not open the input file!");
    }

    // Initalize the variables for the Simulatiom
    int particle_number;
    double time_step;
    double total_time;

    // Initalize the variables for the Diagnostic
    std::string *params = new std::string[9];
    double *binsize = new double[9];
    double *binmax = new double[9];
    double *binmin = new double[9];
    int n_par = 0;
    double* bin_n = new double[9];

    // Setup the variables
    setupInputVariable(input_file, particle_number, time_step, total_time, params, binsize, binmax, binmin, bin_n, n_par);
   
    // ============================
    //   START OF MAIN CODE BLOCK
    // ============================

    Particle* particles = new Particle[particle_number];  // Array of pointers

    createParticles(particles, particle_number);

    for (double t = 0; t <= total_time; t += time_step){

        std::vector<int>* histograms = new std::vector<int>[n_par]; // For Diagnostics Purposes

        for (int i = 0; i < n_par; i++){
            std::cout <<"a  "<<bin_n[i]<<std::endl;
            histograms[i] = std::vector<int>(bin_n[i], 0); // Create the histograms necessary
        }
        
        // Testing the Diagnostics

        PerformDiagnostics(histograms, particles[0], params, binsize, binmax, binmin, n_par); // use in the boris pusher

        for (int value : histograms[0]) {
            std::cout << value << " ";
        }

    }

	// ===========================
	//   PRINTS TO CHECK HEALTH
	// ===========================

    //particles[0].display_position();
    //std::cout<< particles[0].getPosition()[1] <<std::endl;
    //double newpos[3] = {1,1,2};
    //double teste[3] = {1,0,1};

    //double* a = cross(newpos, teste);
    //double b = inner(newpos, teste);
    //std::cout << "Cross product: "<<a[0]<<", "<<a[1]<<", "<<a[2]<<std::endl;
    //std::cout << "Inner Product: "<<b<<std::endl;
    
    //particles[0].setPosition(newpos);
    //std::cout<< particles[0].getPosition()[1] <<std::endl;

    // testing omega
    double u[3] = {1,0,0};
    double E[3] = {0,1,0};
    double B[3] = {0,0,1};

    double* om = Omega(1.41, u , E, B);
    std::cout << "Omega: "<<om[0]<<", "<<om[1]<<", "<<om[2]<<std::endl;
    std::cout <<"lal"<<std::endl;
    

	return 0;
}