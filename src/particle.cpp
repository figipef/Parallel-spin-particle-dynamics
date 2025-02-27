#include "particle.hpp"

// Default Consctructor
Particle::Particle(){}

// Constructor
Particle::Particle(double pos[3], double mom[3], double spn[3], int tg) : tag(tg) {
    std::copy(pos, pos + 3, position);
    std::copy(mom, mom + 3, momentum);
    std::copy(spn, spn + 3, spin);
}
// Member functions to set values
void Particle::setPosition(double pos[3]){
    std::copy(pos, pos + 3, position);
}
void Particle::setMomentum(double mom[3]){
    std::copy(mom, mom + 3, momentum);
}
void Particle::setSpin(double spn[3]){
    std::copy(spn, spn + 3, spin);

}
// Member functions to get values
const double* Particle::getPosition() const{
    return position;
}
const double* Particle::getMomentum() const{
    return momentum;
}
const double* Particle::getSpin() const{
    return spin;
}
const int Particle::getTag() const{
    return tag;
}
// Functions to display values
void Particle::display_position() const{
    std::cout << "[" << position[0] << ", " << position[1] << ", " << position[2] << "]\n";
}
void Particle::display_momentum() const{
    std::cout << "[" << momentum[0] << ", " << momentum[1] << ", " << momentum[2] << "]\n";
}
void Particle::display_spin() const{
    std::cout << "[" << spin[0] << ", " << spin[1] << ", " << spin[2] << "]\n";
}
