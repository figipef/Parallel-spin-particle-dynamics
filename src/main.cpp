#include "particle.hpp"
#include "functions.hpp"
#include "laser.hpp"

#include <random>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <vector>

int main() {

    // ============================
    //   VARIABLE INTIALIZATION
    // ============================

	std::ifstream input_file("input.txt"); // Open the input file

    if (!input_file) {
        throw std::runtime_error("Could not open the input file!");
    }

    // Initalize the variables for the Simulatiom
    
    double time_step;
    double total_time;

    // Initalize the variables for particle creation

    int particle_number;
    std::string *distribution_types = new std::string[3]; // "uni" or "gau"
    double *distribution_sizes = new double[3]; // size of uniforms and std for gaussians

    // Initalize the variables for the Diagnostic
    std::string *params = new std::string[9];
    double *binsize = new double[9];
    double *binmax = new double[9];
    double *binmin = new double[9];
    double *fieldiag = new double[6];
    double *mom_dir = new double[3];
    double *spin_dir = new double[3];
    int RR = 0;
    int n_par = 0;
    int* bin_n = new int[9];

    int step_diag;

    int* follow_params = new int[3];

    // Initaliaze the variables for the lasers

    Laser* lasers = new Laser[10]; // Not the real amount of Lasers, just allocating memory for future use
    int laser_number;

    // Setup the variables
    setupInputVariable(input_file, particle_number, distribution_types, distribution_sizes, mom_dir, spin_dir, time_step, total_time, step_diag, params, binsize, binmax, binmin, bin_n, n_par, fieldiag, lasers, laser_number, RR, follow_params);

    DiagnosticParameters diag_params(params, binsize, binmax, binmin, n_par); // Save the diagnostics to a struct for easier usage

    //std::ofstream file_pos("../output/position.txt");
    //std::ofstream file_mom("../output/momentum.txt");
    //std::ofstream file_spn("../output/spin.txt");

    std::ofstream file_electric_y("../output/e_field.txt");

    // ============================
    //   START OF MAIN CODE BLOCK
    // ============================

    // Create the particle array
    Particle* particles = new Particle[particle_number];  // Array of pointers

    createParticles(particles, particle_number, distribution_types, distribution_sizes, mom_dir, spin_dir, lasers, laser_number);

    int* followed_particles;
    int number_follow_particles = 0;

    std::ofstream file_pos[100];
    std::ofstream file_mom[100];
    std::ofstream file_spn[100];

    if (follow_params[2] > particle_number){std::cout <<" Number of Particles is too high\n";}

    if (follow_params[0] == 1){

        if (follow_params[1] == 1){

            followed_particles = new int[follow_params[2]];

            std::random_device rd;
            std::mt19937 gen(rd()); // Mersenne Twister engine
            std::uniform_int_distribution<> dis(0, particle_number);

            number_follow_particles = follow_params[2];

            for (int i = 0; i < number_follow_particles; i++){

                followed_particles[i] = dis(gen);
                std::cout << followed_particles[i] << "\n";
                
                file_pos[i].open("../output/position" + std::to_string(followed_particles[i]) + ".txt");
                file_mom[i].open("../output/momentum" + std::to_string(followed_particles[i]) + ".txt");
                file_spn[i].open("../output/spin" + std::to_string(followed_particles[i]) + ".txt");

            }

        } else if (follow_params[1] == 0){

            followed_particles = new int[1];

            number_follow_particles = 1;

            followed_particles[0] = follow_params[2];

            file_pos[0].open("../output/position" + std::to_string(followed_particles[0]) + ".txt");
            file_mom[0].open("../output/momentum" + std::to_string(followed_particles[0]) + ".txt");
            file_spn[0].open("../output/spin" + std::to_string(followed_particles[0]) + ".txt");

        } else {

            std::cout <<"smth is wrong with following \n";
        }
    }

    int counter = 0; // Counter for the diagnostics file numbering

    for (double t = 0; t <= total_time; t += time_step){

        printProgressBar(t, total_time, 50);

        // Perform Diagnostics

        if (step_diag >= 1 && counter % step_diag == 0  && n_par > 0){

            Histogram hist = createHistogram(n_par, bin_n);

            boris(particles, lasers, t, time_step, particle_number, laser_number, RR, &hist, &diag_params);
            
            writeDiagnosticsToFile(hist, counter, t, params, n_par);

        } else { // Normal Boris run

            boris(particles, lasers, t, time_step, particle_number, laser_number, RR);
        }

        if (follow_params[0] == 1){

            for (int i = 0; i < number_follow_particles; i++){
                writeToFile(file_pos[i], particles[followed_particles[i]], 'p');
                writeToFile(file_mom[i], particles[followed_particles[i]], 'm');
                writeToFile(file_spn[i], particles[followed_particles[i]], 's');
            }                
        }

        // Following the information on particle 0
        //writeToFile(file_pos, particles[0], 'p');
        //writeToFile(file_mom, particles[0], 'm');
        //writeToFile(file_spn, particles[0], 's');

        counter++;
    }

	// ===========================
	//   PRINTS TO CHECK HEALTH
	// ===========================
    /*
    // Lasers health
    std::cout << "---Lasers--- \n";
    for (int i = 0; i < laser_number; i++){
        std::cout <<"Laser number = " << i <<std::endl;
        std::cout << "Tag = " << lasers[i].tag << "\n";
        std::cout << "Type = " << lasers[i].type << "\n";
        std::cout << "Freq = " << lasers[i].freq << "\n";
        std::cout << "env_freq = " << lasers[i].env_freq << "\n";
        std::cout << "length = " << lasers[i].length << "\n";
        std::cout << "E = [ " <<lasers[i].get_E_0()[0]<<", "<<lasers[i].get_E_0()[1]<<", "<<lasers[i].get_E_0()[2]<< " ]\n";
        std::cout << "B = [ " <<lasers[i].get_B_0()[0]<<", "<<lasers[i].get_B_0()[1]<<", "<<lasers[i].get_B_0()[2]<< " ]\n";
    }

    std::cout << "----Testing----\n";
    double t = 0;
    int tipo = 3;
    double pos[3] = {1,0,0};

    for (int i = 0; i<6;i++){
        
        std::cout << " Component "<<i+1<<" : "<< lasers[0].get_fields(&t, pos)[i]<<"\n";
    }

    //std::cout <<"Resultado: "<<cos(pos[0]*lasers[0].k[0] + pos[1]*lasers[0].k[1]+pos[2]*lasers[0].k[2] - t * 1.41421) * sin((pos[0]*lasers[0].k[0] + pos[1]*lasers[0].k[1]+pos[2]*lasers[0].k[2] - t * 2)/5) * sin((pos[0]*lasers[0].k[0] + pos[1]*lasers[0].k[1]+pos[2]*lasers[0].k[2] - t * 2)/5)<<"\n";
    std::cout << "---------\n\n";
    
    std::cout<< "Updating, 1 step" <<std::endl;

    particles[0].display_position();
    particles[0].display_momentum();
    particles[0].display_spin();

    std::cout << " Field diagnostics test:\nShow Efield " << *fieldiag << "\nShow Bfield " << fieldiag[1] << "\nBin size " << fieldiag[2] << "\nBin start " << fieldiag[3] << "\nBin end " << fieldiag[4] << "\n";

    double dt = 1;
    int iter = 30;
    FieldDiagWritter(dt, iter, fieldiag, lasers, laser_number);

    int N = particle_number;
    std::ofstream file_lots_spin("../output/lots_spin.txt");
    for (int n = 0; n < N; n++){
        writeToFile(file_lots_spin, particles[n], 's');
    }
    */
    std::cout <<"\nshould be finished\n";

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
    delete[] mom_dir;
    delete[] spin_dir;

    delete[] follow_params;
    delete[] followed_particles;

	return 0;
}