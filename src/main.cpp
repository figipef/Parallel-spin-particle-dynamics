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

	std::ifstream input_file("input.txt"); // Open the input file

    if (!input_file) {
        std::cerr << "Error: Could not open the file!" << std::endl;
        return 1;
    }

    std::unordered_map<std::string, std::string> values;
    std::string line;
    
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

    int particle_number = std::stoi(values["NUMBER_OF_PARTICLES"]);
    double time_step = std::stod(values["TIME_STEP"]);
    int total_time = std::stoi(values["TOTAL_TIME"]);

    std::string *params = new std::string[3];
    double *binsize = new double[3];
    double *binmax = new double[3];
    double *binmin = new double[3];
    int n_par = 0;

    // Save the Diagnostics parameters to strings
    params[0] = values["PAR1"];
    params[1] = values["PAR2"];
    params[2] = values["PAR3"];

    // Save the wanted Bin data
    
    // Save the BinSize values to an array
    if (!values["BIN_SIZE_1"].empty()){
        binsize[0] = std::stod(values["BIN_SIZE_1"]);
        n_par += 1;
    }
    if (!values["BIN_SIZE_2"].empty()){
        binsize[1] = std::stod(values["BIN_SIZE_2"]);
        n_par += 1;
    }
    if (!values["BIN_SIZE_3"].empty()){
        binsize[2] = std::stod(values["BIN_SIZE_3"]);
        n_par += 1;
    }

    // Save the BinMax values to an array
    if (!values["B1_MAX"].empty()){
        binmax[0] = std::stod(values["B1_MAX"]);
    }
    if (!values["B2_MAX"].empty()){
        binmax[1] = std::stod(values["B2_MAX"]);
    }
    if (!values["B3_MAX"].empty()){
        binmax[2] = std::stod(values["B3_MAX"]);
    }

    // Save the BinMin values to an array 
    if (!values["B1_MIN"].empty()){
        binmin[0] = std::stod(values["B1_MIN"]);
    }
    if (!values["B2_MIN"].empty()){
        binmin[1] = std::stod(values["B2_MIN"]);
    }
    if (!values["B3_MIN"].empty()){
        binmin[2] = std::stod(values["B3_MIN"]);
    }

    double* bin_n = new double[n_par]; 
    for (int i = 0; i < n_par; i++){
        bin_n[i] = (binmax[i] - binmin[i]) / binsize[i] + 1; 
    }

    // ============================
    //   START OF MAIN CODE BLOCK
    // ============================

    Particle* particles = new Particle[particle_number];  // Array of pointers

    createParticles(particles, particle_number);

    for (double t = 0; t <= total_time; t += time_step){

        std::vector<int>* histograms = new std::vector<int>[n_par]; // For Diagnostics Purposes

        for (int i = 0; i < n_par; i++){
            histograms[i] = std::vector<int>(bin_n[i], 0); // Create the histograms necessary
        }

        // Testing the Diagnostics


        PerformDiagnostics(histograms, particles[0], params, binsize, binmax, binmin, n_par);

        for (int value : histograms[0]) {
            std::cout << value << " ";
        }

    }

	// ===========================
	//   PRINTS TO CHECK HEALTH
	// ===========================

    particles[0].display_position();
    std::cout<< particles[0].getPosition()[1] <<std::endl;
    double newpos[3] = {0,0.5,0};
    particles[0].setPosition(newpos);
    std::cout<< particles[0].getPosition()[1] <<std::endl;

	return 0;
}