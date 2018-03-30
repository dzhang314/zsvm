#ifndef ZSVM_SPHERICAL_ECG_PAIRWISE_CONTEXT_HPP_INCLUDED
#define ZSVM_SPHERICAL_ECG_PAIRWISE_CONTEXT_HPP_INCLUDED

#include <cmath>
#include <string>
#include <vector>

#include "EnumTypes.hpp"
#include "Particle.hpp"
#include "PackedLinearAlgebra.hpp"

namespace zsvm {

    class SphericalECGPairwiseContext {

        const std::size_t num_particles;
        const std::vector<double> adjusted_carrier_values;

        static std::vector<double> extract_adjusted_carrier_values(
                long long int space_dimension,
                const std::vector<Particle> &particles,
                const PairwisePotential &pairwise_potential) {
            std::vector<double> result;
            const std::string &carrier_name = pairwise_potential.carrier;
            const std::size_t n = particles.size();
            const auto d = static_cast<double>(space_dimension);
            const double dimension_factor =
                    std::tgamma(0.5 * (d - 1.0)) / std::tgamma(0.5 * d);
            for (std::size_t i = 0; i < n - 1; ++i) {
                for (std::size_t j = i + 1; j < n; ++j) {
                    const auto carrier_i_iter =
                            particles[i].carriers.find(carrier_name);
                    if (carrier_i_iter == particles[i].carriers.end()) {
                        std::cerr << "ERROR: Particle lacks required carrier "
                                  << carrier_name << "." << std::endl;
                        std::exit(EXIT_FAILURE);
                    }
                    const auto carrier_j_iter =
                            particles[j].carriers.find(carrier_name);
                    if (carrier_j_iter == particles[j].carriers.end()) {
                        std::cerr << "ERROR: Particle lacks required carrier "
                                  << carrier_name << "." << std::endl;
                        std::exit(EXIT_FAILURE);
                    }
                    result.push_back(dimension_factor
                                     * pairwise_potential.strength
                                     * carrier_i_iter->second
                                     * carrier_j_iter->second);
                }
            }
            return result;
        }

        explicit SphericalECGPairwiseContext(
                long long int space_dimension,
                const std::vector<Particle> &particles,
                const PairwisePotential &pairwise_potential)
                : num_particles(particles.size()),
                  adjusted_carrier_values(extract_adjusted_carrier_values(
                          space_dimension, particles, pairwise_potential)) {}

        double evaluate_kernel(const double *__restrict__ d) {
            double result = 0.0;
            for (std::size_t i = 0, k = 0; i < num_particles - 1; ++i) {
                for (std::size_t j = i + 1; j < num_particles; ++j, ++k) {
                    const double alpha = d[i * (i + 3) / 2] + d[j * (j + 3) / 2]
                                         - 2.0 * d[j * (j + 1) / 2 + i];
                    result += adjusted_carrier_values[k] /
                              std::sqrt(2.0 * alpha);
                }
            }
            return result;
        }

    }; // class SphericalECGPairwiseContext

} // namespace zsvm

#endif // ZSVM_SPHERICAL_ECG_PAIRWISE_CONTEXT_HPP_INCLUDED
