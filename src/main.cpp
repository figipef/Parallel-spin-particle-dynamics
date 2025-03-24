#include "particle.hpp"
#include "functions.hpp"
#include "laser.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <vector>

#include <mpi.h>

int main(int argc, char** argv) {

    // ============================
    //   VARIABLE INTIALIZATION
    // ============================

	std::ifstream input_file("input.txt"); // Open the input file

    if (!input_file) {
        std::cerr << "Could not open the input file!" << std::endl;
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
    double *spin_dir = new double[3];
    int n_par = 0;
    int* bin_n = new int[9];

    int step_diag;

    // Initaliaze the variables for the lasers

    Laser* lasers = new Laser[10]; // Not the real amount of Lasers, just allocating memory for future use
    int laser_number;

    // Setup the variables
    setupInputVariable(input_file, particle_number, distribution_types, distribution_sizes, spin_dir, time_step, total_time, step_diag, params, binsize, binmax, binmin, bin_n, n_par, fieldiag, lasers, laser_number);

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

    createParticles(particles, particle_number, distribution_types, distribution_sizes, spin_dir, lasers, laser_number);

    // Divide the particle number to each process

    // ==================================
    //    PARALELIZATION INITIALIZATION
    // ==================================

    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double start_time = MPI_Wtime(); // Performance Time check

    int local_particle_count = particle_number / size;
    int start_idx = rank * local_particle_count;
    int end_idx = (rank == size - 1) ? particle_number : start_idx + local_particle_count;

    // Prints for health
    /*
    std::cout << "Electric field 0: "<<lasers[0].get_E_0()[0]<<", "<<lasers[0].get_E_0()[1]<<", "<<lasers[0].get_E_0()[2]<<std::endl;

    std::cout << "Magnetic field 0: "<<lasers[0].get_B_0()[0]<<", "<<lasers[0].get_B_0()[1]<<", "<<lasers[0].get_B_0()[2]<<std::endl;

    particles[0].display_position();
    particles[0].display_momentum();
    particles[0].display_spin();
    */
    int counter = 0; // Counter for the diagnostics file numbering
    for (double t = 0; t <= total_time; t += time_step){

        // Perform Diagnostics

        if (step_diag >= 1 && counter % step_diag == 0 && n_par > 0){

            Histogram local_hist = createHistogram(n_par, bin_n);
            Histogram global_hist = createHistogram(n_par, bin_n);

            boris(particles + start_idx, lasers, t, time_step, end_idx - start_idx, laser_number, &local_hist, &diag_params);

            int total_bins = (local_hist.is_matrix) ? (bin_n[0] * bin_n[1]) : bin_n[0];

            // Flatten the histogram for communication
            std::vector<int> local_flat = local_hist.flatten();
            std::vector<int> global_flat(total_bins, 0);

            // Perform MPI_Reduce
            MPI_Reduce(local_flat.data(), global_flat.data(), total_bins, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

            // Rank 0 reconstructs global histogram
            if (rank == 0) {
                global_hist.from_flattened(global_flat, bin_n);
                writeDiagnosticsToFile(global_hist, counter, t);
            }

        } else { // Normal Boris run

            boris(particles + start_idx, lasers, t, time_step, end_idx - start_idx, laser_number);
        }

        // Following the information on particle 0
        if (rank == 0) {
            writeToFile(file_pos, particles[0], 'p');
            writeToFile(file_mom, particles[0], 'm');
            writeToFile(file_spn, particles[0], 's');
        }
        counter++;
    }


    // Calculate the time taken to process everything

    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;
    double max_time;
    MPI_Reduce(&elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << "Max Execution Time across all processes: " << max_time << " seconds\n";
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

    int N = 10000;
    std::ofstream file_lots_spin("../output/lots_spin.txt");
    for (int n = 0; n < N; n++){
        writeToFile(file_lots_spin, particles[n], 's');
    }

    std::cout <<"should be finished\n";

    */
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
    delete[] spin_dir;

    MPI_Finalize();
	return 0;
}