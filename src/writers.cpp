#include "writers.hpp"

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
 
void writeDiagnosticsToFile(const Histogram hist, const int iter, const double t, const std::string* params, const int n_of_par){

	std::ostringstream filename_temp;

	filename_temp << "../output/histogram" << std::setw(10) << std::setfill('0') << iter;  // 8 digits with leading zeros

	std::string filename = filename_temp.str() + ".txt";

	std::ofstream outFile(filename);

	if (outFile.is_open()){
		outFile << "Time: " << t;

		for (int i; i < n_of_par; i++){
			outFile << " Par" << i+1 <<": "<< params[i];
		}
		outFile << "\n";
		if (!hist.is_matrix){

			 for (const auto& bin : hist.vector1D) {
        	    outFile << bin << "\n";  // Each bin value on a new line
        	}
        	outFile.close();
        	std::cout << "1D Histogram saved to " << filename << std::endl;

		} else {

			for (const auto& row : hist.matrix2D) {

        	    for (size_t i = 0; i < row.size(); ++i) {

        	        outFile << row[i];
        	        if (i < row.size() - 1) outFile << "\t";  // Tab delimiter between columns

        	    }

        	    outFile << "\n";  // Newline for each row
        	}

        	outFile.close();
        	std::cout << "2D Histogram saved to " << filename << std::endl;
		}

	} else {
		std::cerr << "Unable to open the file for writing!" << std::endl;
	}
}

// Helper function to help parse vectors to strings
void parseVector(const std::string& value, double vec[3]) {
    std::istringstream ss(value);
    char comma; // To skip the commas
    ss >> vec[0] >> comma >> vec[1] >> comma >> vec[2];
}

// Courtesy of João Chaveiro
void printProgressBar(int progress, int total, int barWidth) {

  float percent = (float)progress / total;
  int filled = percent * barWidth;

  std::cout << "\r[";
  for (int i = 0; i <= barWidth; i++) {
      if (i < filled)
          std::cout << "=";
      else if (i == filled)
          std::cout << ">";
      else
          std::cout << " ";
  }
  std::cout << "] " << int(percent * 100.0) << "%";
  std::cout.flush();
}