#ifndef ZSVM_SPHERICAL_ECG_CONFINING_CONTEXT_HPP_INCLUDED
#define ZSVM_SPHERICAL_ECG_CONFINING_CONTEXT_HPP_INCLUDED

#include <string>
#include <vector>

#include "EnumTypes.hpp"
#include "Particle.hpp"

namespace zsvm {

    class SphericalECGConfiningContext { // TODO

    public:

        explicit SphericalECGConfiningContext(
                const std::vector<zsvm::Particle<double>> &,
                const ConfiningPotential<double> &) {}

        double evaluate_kernel() const noexcept {
            return 0.0;
        }

    }; // class SphericalECGConfiningContext

} // namespace zsvm

#endif // ZSVM_SPHERICAL_ECG_CONFINING_CONTEXT_HPP_INCLUDED
