#include "particle.hpp"
#include "functions.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>

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

    // Save the Diagnostics parameters to strings
    std::string par1 = values["PAR1"];
    std::string par2 = values["PAR2"];
    std::string par3 = values["PAR3"];

    // Save the wanted Bin data
    double bin1_size = std::stod(values["BIN_SIZE_1"]);
    double bin2_size = std::stod(values["BIN_SIZE_2"]);
    //double bin3_size = std::stod(values["BIN_SIZE_3"]);

    double bin1_max = std::stod(values["B1_MAX"]);
    double bin2_max = std::stod(values["B2_MAX"]);
    //double bin3_max = std::stod(values["B3_MAX"]);

    double bin1_min = std::stod(values["B1_MIN"]);
    double bin2_min = std::stod(values["B2_MIN"]);
    //double bin3_min = std::stod(values["B3_MIN"]);


    // ============================
    //   START OF MAIN CODE BLOCK
    // ============================

    Particle* particles = new Particle[particle_number];  // Array of pointers
    /*
    double mom[3] = {0,0,0};
    double pos[3] = {0,0,0};
    double spin[3] = {0,1,0};

	for (int i = 0; i < particle_number; ++i) {
	    particles[i] = new Particle(mom,pos,spin,i);  // Dynamically allocate objects
	}
    */
    createParticles(particles, particle_number);
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