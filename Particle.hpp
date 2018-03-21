#ifndef ZSVM_PARTICLE_HPP_INCLUDED
#define ZSVM_PARTICLE_HPP_INCLUDED

// C++ standard library headers
#include <iostream> // for std::ostream
#include <map> // for std::map
#include <string> // for std::string

namespace zsvm {

    struct Particle {

        const std::size_t type_id;
        const int spin;
        const std::map<std::string, double> carriers;

        Particle(std::size_t type_id, int spin,
                 const std::map<std::string, double> &carriers);

    }; // struct Particle

} // namespace zsvm

#endif // ZSVM_PARTICLE_HPP_INCLUDED
