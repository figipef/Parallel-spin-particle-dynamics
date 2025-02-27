#pragma once

class Particle {
private:

    double position[3];
    double velocity[3];
    double spin[3];
    int tag;

public:
    // Default Constructor
    Particle();
    // Constructor
    Particle(double[3], double[3], double[3], int);

    // Member functions to set values
    void setPosition(double[3]);
    void setVelocity(double[3]);
    void setSpin(double[3]);

    // Member functions to get values
    double* getPosition() const;
    double* getVelocity() const;
    double* getSpin() const;
    int getTag() const;

    // Functions to display values
    void display_positon() const;
    void display_velocity() const;
    void display_spin() const;

};
