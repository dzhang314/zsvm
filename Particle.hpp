#ifndef ZSVM_PARTICLE_HPP_INCLUDED
#define ZSVM_PARTICLE_HPP_INCLUDED

// C++ standard library headers
#include <map>
#include <utility> // for std::move
#include <string>

namespace zsvm {

    template <typename T>
    struct Particle {

        const std::size_t type_id;
        const int spin;
        const std::map<std::string, T> carriers;

        Particle(std::size_t type_id, int spin,
                 std::map<std::string, T> carriers) :
                type_id(type_id), spin(spin), carriers(std::move(carriers)) {}

    }; // struct Particle

} // namespace zsvm

#endif // ZSVM_PARTICLE_HPP_INCLUDED
