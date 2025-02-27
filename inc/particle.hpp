#pragma once

#include <algorithm>
#include <iostream>

class Particle {
private:

    double position[3];
    double momentum[3];
    double spin[3];
    int tag;

public:
    // Default Constructor
    Particle();
    // Constructor
    Particle(double[3], double[3], double[3], int);

    // Member functions to set values
    void setPosition(double[3]);
    void setMomentum(double[3]);
    void setSpin(double[3]);

    // Member functions to get values
    const double* getPosition() const;
    const double* getMomentum() const;
    const double* getSpin() const;
    const int getTag() const;

    // Functions to display values
    void display_position() const;
    void display_momentum() const;
    void display_spin() const;

};
