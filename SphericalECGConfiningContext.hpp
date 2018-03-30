#ifndef ZSVM_SPHERICAL_ECG_CONFINING_CONTEXT_HPP_INCLUDED
#define ZSVM_SPHERICAL_ECG_CONFINING_CONTEXT_HPP_INCLUDED

#include <string>
#include <vector>

#include "EnumTypes.hpp"
#include "Particle.hpp"

namespace zsvm {

    class SphericalECGConfiningContext { // TODO

        explicit SphericalECGConfiningContext(
                const std::vector<Particle> &particles,
                const ConfiningPotential &confining_potential) {}

        double evaluate_kernel() {
            return 0.0;
        }

    }; // class SphericalECGConfiningContext

} // namespace zsvm

#endif // ZSVM_SPHERICAL_ECG_CONFINING_CONTEXT_HPP_INCLUDED
