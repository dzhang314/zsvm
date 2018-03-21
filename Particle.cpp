#include "Particle.hpp"


zsvm::Particle::Particle(std::size_t type_id, int spin,
                         const std::map<std::string, double> &carriers)
        : type_id(type_id), spin(spin), carriers(carriers) {}
