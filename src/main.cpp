#include "particle.hpp"
#include "functions.hpp"
#include "laser.hpp"
#include "logger.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <chrono>

int main(int argc, char* argv[]) {

    // ============================
    //   VERBOSITY INTIALIZATION
    // ============================

    int verbosity = 0; // Default verbosity level

    if (argc > 1) {
        verbosity = std::atoi(argv[1]); // Convert argument to integer
    }

    Logger logger(verbosity);

    logger.log(0,"\n[0] Starting Program \n");

    // ============================
    //   VARIABLE INTIALIZATION
    // ============================


    // Initalize the variables for the Simulatiom
    
    double time_step;                // Simulation time-step
    double total_time;               // Total time the simulation is ran for
    int RR = 0;                      // Bool to check for the inclusion of radiation reaction
    int* follow_params = new int[3]; // Refer to read me to see what each index mean

    // Initalize the variables for particle creation

    int particle_number;
    std::string *distribution_types = new std::string[3]; // Refer to read me to see the type and stored values
    double *distribution_sizes = new double[3];           // Refer to read me to see the size values
    double *pos_dir = new double[3];                      // Prefered inital position
    double *mom_dir = new double[3];                      // Prefered momentum direction
    double *spin_dir = new double[3];                     // Prefered spin direction

    // Initalize the variables for the Diagnostic

    std::string *params = new std::string[9]; // All parameters in the form "p1" or "m3" etc
    double *binsize = new double[9];          // The size of bin for each parameter
    double *binmax = new double[9];           // Largest values recorded for each parameter
    double *binmin = new double[9];           // Lowest values recorded for each parameter
    int* bin_n = new int[9];                  // Number of bins for each parameter
    double *fieldiag = new double[6];         // Parameters for field diagnostics
    int n_par = 0;                            // Number of parameters
    int step_diag;                            // Number of iterations for diagnostics

    // Initialize the variables to follow particles

    int* followed_particles;          // Array with the indexes for the followed particless
    int number_follow_particles = 0;  // Number of followed particles

    std::ofstream file_pos[100];      // Array to the files to store the position of the particles
    std::ofstream file_mom[100];      // Array to the files to store the momentum of the particles 
    std::ofstream file_spn[100];      // Array to the files to store the spin of the particles  

    // Initaliaze the variables for the lasers

    Laser* lasers = new Laser[10]; // Not the real amount of Lasers, just allocating memory for future use
    int laser_number;              // Number of lasers for easir looping

    // ============================
    //       SIMULATION SETUP
    // ============================

    // Start the setup timer
    auto time_setup = std::chrono::high_resolution_clock::now();

    logger.log(1,"\n [1] Reading Input File \n");

    std::ifstream input_file("input.txt"); // Open the input file

    if (!input_file) { throw std::runtime_error("Input file not found!"); } // Check the existence of the input file

    // Setup ALL the input file variables
    setupInputVariable(input_file, particle_number, distribution_types, distribution_sizes, pos_dir, mom_dir, spin_dir, time_step, total_time, step_diag, params, binsize, binmax, binmin, bin_n, n_par, fieldiag, lasers, laser_number, RR, follow_params);

    logger.logLasers(3, lasers, laser_number);

    // Save the diagnostics to a struct for easier usage 
    DiagnosticParameters diag_params(params, binsize, binmax, binmin, n_par); 

    logger.logDiag(3, diag_params); 

    Particle* particles = new Particle[particle_number];  // Create the particle array

    logger.log(1,"\n [1] Creating Particle Array \n");    

    // Create the particles according to input parameters
    createParticles(particles, particle_number, distribution_types, distribution_sizes, pos_dir, mom_dir, spin_dir, lasers, laser_number);
    
    logger.log(2,"\n  [2] Setting up followed particles \n"); 

    // Create the files and save the necessary variables to follow particles
    setupFollowParticles(follow_params, file_pos, file_mom, file_spn, followed_particles, number_follow_particles, particle_number);

    // ============================
    //     MAIN LOOP CODE BLOCK
    // ============================

    // Start the Loop timer
    auto time_loop = std::chrono::high_resolution_clock::now();

    logger.log(1,"\n [1] Starting Simulation \n");    

    int counter = 0; // Counter for the diagnostics file numbering

    for (double t = 0; t <= total_time; t += time_step){

        printProgressBar(t, total_time, 50); // Dinamic progress bar printing

        if (step_diag >= 1 && counter % step_diag == 0  && n_par > 0){

            // Perform Diagnostics

            Histogram hist = createHistogram(n_par, bin_n);

            boris(particles, lasers, t, time_step, particle_number, laser_number, RR, &hist, &diag_params);
            
            writeDiagnosticsToFile(hist, counter, t, params, n_par);

        } else { // Normal Boris run

            boris(particles, lasers, t, time_step, particle_number, laser_number, RR);
        }

        if (follow_params[0] == 1){

            // Save write to the files the necessary particle information

            for (int i = 0; i < number_follow_particles; i++){
                writeToFile(file_pos[i], particles[followed_particles[i]], 'p');
                writeToFile(file_mom[i], particles[followed_particles[i]], 'm');
                writeToFile(file_spn[i], particles[followed_particles[i]], 's');
            }                
        }

        counter++;
    }
    std::cout <<"\n";

    // Finish the Loop timer
    auto time_loop_end = std::chrono::high_resolution_clock::now();

    logger.log(1,"\n [1] Successfully Finished Simulation \n");

    logger.log(2,"\n  [2] Writing Fields Diagnostics (If selected) \n");
    FieldDiagWritter(time_step, counter, fieldiag, lasers);

    //  Clean the Dinamically allocated memory

    delete[] lasers;
    delete[] particles;

    delete[] distribution_types;
    delete[] distribution_sizes;
    delete[] params;
    delete[] binsize;
    delete[] binmax;
    delete[] binmin;
    delete[] bin_n;
    delete[] fieldiag;
    delete[] pos_dir;
    delete[] mom_dir;
    delete[] spin_dir;

    delete[] follow_params;
    delete[] followed_particles;

    auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(time_loop - time_setup);
    auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(time_loop_end - time_loop);

    logger.log(2,"\n  [2] Setup took " + std::to_string(static_cast<double>(duration1.count())/1000.) + " seconds \n");
    logger.log(2,"\n  [2] Simulation took " + std::to_string(static_cast<double>(duration2.count())/1000.) + " seconds \n");

    logger.log(0,"\n[0] Successfully Finished Running \n");

	return 0;
}