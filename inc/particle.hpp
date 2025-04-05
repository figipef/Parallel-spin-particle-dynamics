#pragma once

#include <algorithm>
#include <iostream>

class Particle {

private:

    double position[3]; // 3D position (x,y,z)
    double momentum[3]; // 3D momentum (px,py,pz)
    double spin[3];     // 3D spin (Sx,Sy,Sz)
    int tag;            // Particle Tag

public:

    // Default Constructor
    Particle();
    // Constructor
    Particle(double[3], double[3], double[3], int); // initial position, momento, spin and tag

    // Member functions to set values
    void setPosition(double[3]);
    void setMomentum(double[3]);
    void setSpin(double[3]);

    // Member functions to get values
    const double* getPosition() const;
    const double* getMomentum() const;
    const double* getSpin() const;
    int getTag() const;

    // Functions to display values
    void display_position() const;
    void display_momentum() const;
    void display_spin() const;

};
