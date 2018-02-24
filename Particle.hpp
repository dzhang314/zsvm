#ifndef ZSVM_PARTICLE_HPP_INCLUDED
#define ZSVM_PARTICLE_HPP_INCLUDED

// C++ standard library headers
#include <iostream>

namespace zsvm {

    enum class Spin {
        UP,
        DOWN,
    };

    struct Particle {
        int type;
        double mass;
        double charge;
        Spin spin;
    };

} // namespace zsvm

std::ostream &operator<<(std::ostream &os, zsvm::Spin spin);

std::ostream &operator<<(std::ostream &os, const zsvm::Particle &particle);

#endif // ZSVM_PARTICLE_HPP_INCLUDED
