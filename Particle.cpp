#include "Particle.hpp"

std::ostream &operator<<(std::ostream &os, zsvm::Spin spin) {
    switch (spin) {
        case zsvm::Spin::UP:
            os << "UP";
            break;
        case zsvm::Spin::DOWN:
            os << "DOWN";
            break;
    }
    return os;
}

std::ostream &operator<<(std::ostream &os, const zsvm::Particle &particle) {
    return os << "Particle(type=" << particle.type
              << ", mass=" << particle.mass
              << ", charge=" << particle.charge
              << ", spin=" << particle.spin << ")";
}
