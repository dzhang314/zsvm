#ifndef ZSVM_SPHERICAL_ECG_DISPERSION_CONTEXT_HPP_INCLUDED
#define ZSVM_SPHERICAL_ECG_DISPERSION_CONTEXT_HPP_INCLUDED

#include <string>
#include <vector>

#include "EnumTypes.hpp"
#include "Particle.hpp"
#include "PackedLinearAlgebra.hpp"

namespace zsvm {

    class SphericalECGDispersionContext {

        const std::vector<double> adjusted_carrier_values;
        const packed_kinetic_trace_function_t packed_kinetic_trace;

        static std::vector<double> extract_adjusted_carrier_values(
                long long int space_dimension,
                const std::vector<Particle> &particles,
                const DispersionRelation &dispersion_relation) {
            std::vector<double> result;
            const std::string &carrier_name = dispersion_relation.carrier;
            for (const auto &particle : particles) {
                const auto carrier_iter = particle.carriers.find(carrier_name);
                if (carrier_iter == particle.carriers.end()) {
                    std::cerr << "ERROR: Particle lacks required carrier "
                              << carrier_name << "." << std::endl;
                    std::exit(EXIT_FAILURE);
                }
                result.push_back(space_dimension
                                 * dispersion_relation.strength
                                 * carrier_iter->second);
            }
            return result;
        }

    public:

        explicit SphericalECGDispersionContext(
                long long int space_dimension,
                const std::vector<Particle> &particles,
                const DispersionRelation &dispersion_relation)
                : adjusted_carrier_values(extract_adjusted_carrier_values(
                space_dimension, particles, dispersion_relation)),
                  packed_kinetic_trace(
                          PACKED_KINETIC_TRACE[particles.size() - 1]) {}

        double evaluate_kernel(const double *__restrict__ a,
                               const double *__restrict__ b,
                               const double *__restrict__ d) const noexcept {
            return packed_kinetic_trace(a, b, d,
                                        adjusted_carrier_values.data());
        }

    }; // class SphericalECGDispersionContext

} // namespace zsvm

#endif // ZSVM_SPHERICAL_ECG_DISPERSION_CONTEXT_HPP_INCLUDED
