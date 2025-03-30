#include <iostream>
#include <string>

#include "laser.hpp"
#include "diagnostics.hpp"

class Logger {
private:
    int verbosity;

public:
    explicit Logger(int level) : verbosity(level) {}

    void log(int level, const std::string& message) const {
        if (level <= verbosity) {
            std::cout << message << std::endl;
        }
    }

    void logLasers(int level, const Laser* lasers, int n_of_lasers) const {
        if (level <= verbosity) {
            std::cout <<"\n   [3] Printing used lasers \n \n" ;
            for (int i = 0; i < n_of_lasers; i++){
                std::cout << "    Laser number: "<< i + 1 <<"\n";
                std::cout << "    Laser type  : "<< lasers[i].type <<"\n";
                std::cout << "    k:  [" << lasers[i].k[0]<<" , "<< lasers[i].k[1]<<" , "<< lasers[i].k[2]<<"]\n";
                std::cout << "    E0: [" << lasers[i].e_0[0]<<" , "<< lasers[i].e_0[1]<<" , "<< lasers[i].e_0[2]<<"]\n";
                std::cout << "    B0: [" << lasers[i].b_0[0]<<" , "<< lasers[i].b_0[1]<<" , "<< lasers[i].b_0[2]<<"]\n";
                if (lasers[i].type >= 1){
                    std::cout << "    Laser frequency  : "<< lasers[i].freq <<"\n";
                }
                if (lasers[i].type >= 3){
                    std::cout << "    Laser length     : "<< lasers[i].length <<"\n";
                }
                
                if (lasers[i].type >= 2 && lasers[i].ext_phase != 0){

                    std::cout << "    Laser extra phase: " << lasers[i].ext_phase<<" π\n";
                }
                std::cout <<"\n";
            }

        }
    }

    void logDiag(int level, const DiagnosticParameters diags){
        if (level <= verbosity) {
            std::cout <<"\n   [3] Printing Parameters \n \n" ;

            for (int i = 0; i < diags.n_of_pars; i++){
                std::cout << "    Parameter number: "<< i + 1 <<"\n";
                std::cout << "    Parameter       : "<< diags.params[i][0] << "\n";
                std::cout << "    Axis            : "<< diags.params[i][1] << "\n";
                std::cout << "    Bin size        : "<< diags.bsize[i] << "\n";
                std::cout << "    Bin maximum     : "<< diags.bmax[i] << "\n";
                std::cout << "    Bin minimum     : "<< diags.bmin[i] << "\n";
                std::cout <<"\n";
            }

        }
    }
};