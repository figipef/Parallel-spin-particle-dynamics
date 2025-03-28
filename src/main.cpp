#include "particle.hpp"
#include "functions.hpp"
#include "laser.hpp"

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
    int RR;
    int n_par = 0;
    int* bin_n = new int[9];

    int step_diag;

    // Initaliaze the variables for the lasers

    Laser* lasers = new Laser[10]; // Not the real amount of Lasers, just allocating memory for future use
    int laser_number;

    // Setup the variables
    setupInputVariable(input_file, particle_number, distribution_types, distribution_sizes, mom_dir, spin_dir, time_step, total_time, step_diag, params, binsize, binmax, binmin, bin_n, n_par, fieldiag, lasers, laser_number, RR);

    DiagnosticParameters diag_params(params, binsize, binmax, binmin, n_par); // Save the diagnostics to a struct for easier usage

    std::ofstream file_pos("../output/position.txt");
    std::ofstream file_mom("../output/momentum.txt");
    std::ofstream file_spn("../output/spin.txt");

    std::ofstream file_electric_y("../output/e_field.txt");

    // ============================
    //   START OF MAIN CODE BLOCK
    // ============================

    // Create the particle array
    Particle* particles = new Particle[particle_number];  // Array of pointers

    createParticles(particles, particle_number, distribution_types, distribution_sizes, mom_dir, spin_dir, lasers, laser_number);

    std::cout << "Electric field 0: "<<lasers[0].get_E_0()[0]<<", "<<lasers[0].get_E_0()[1]<<", "<<lasers[0].get_E_0()[2]<<std::endl;

    std::cout << "Magnetic field 0: "<<lasers[0].get_B_0()[0]<<", "<<lasers[0].get_B_0()[1]<<", "<<lasers[0].get_B_0()[2]<<std::endl;

    particles[0].display_position();
    particles[0].display_momentum();
    particles[0].display_spin();

    int counter = 0; // Counter for the diagnostics file numbering
    for (double t = 0; t <= total_time; t += time_step){

        // Perform Diagnostics

        if (counter % step_diag == 0 && step_diag >= 1 && n_par > 0){

            Histogram hist = createHistogram(n_par, bin_n);

            boris(particles,lasers, t, time_step, particle_number, laser_number, RR, &hist, &diag_params);
            
            writeDiagnosticsToFile(hist, counter, t);

        } else { // Normal Boris run

            boris(particles,lasers, t, time_step, particle_number, laser_number, RR);
        }
        
        for (double x =0; x <=20;x = x + 1){
            double pos[3] = {x,0,0};
            //std::cout << lasers[0].get_fields(pos, &t)[1] <<" ";
            file_electric_y << lasers[0].get_fields(&t,pos)[1] <<" ";
        }
        //std::cout <<"\n";
        file_electric_y << "\n";

        // Following the information on particle 0
        writeToFile(file_pos, particles[0], 'p');
        writeToFile(file_mom, particles[0], 'm');
        writeToFile(file_spn, particles[0], 's');

        counter++;
    }

	// ===========================
	//   PRINTS TO CHECK HEALTH
	// ===========================

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

    int N = 10000;
    std::ofstream file_lots_spin("../output/lots_spin.txt");
    for (int n = 0; n < N; n++){
        writeToFile(file_lots_spin, particles[n], 's');
    }

    std::cout <<"should be finished\n";

	return 0;
}