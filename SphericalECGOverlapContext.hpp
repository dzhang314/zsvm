#ifndef ZSVM_SPHERICAL_ECG_OVERLAP_CONTEXT_HPP_INCLUDED
#define ZSVM_SPHERICAL_ECG_OVERLAP_CONTEXT_HPP_INCLUDED

#include <algorithm> // for std::min, std::max
#include <cstddef> // for std::size_t
#include <map> // for std::map
#include <string> // for std::string
#include <vector> // for std::vector

#include "Particle.hpp"
#include "Permutation.hpp"
#include "SphericalECGDispersionContext.hpp"
#include "SphericalECGConfiningContext.hpp"
#include "SphericalECGPairwiseContext.hpp"

namespace zsvm {

    class SphericalECGOverlapContext {

        static std::vector<SphericalECGDispersionContext>
        extract_dispersion_contexts(
                long long int space_dimension,
                const std::vector<Particle<double>> &particles,
                const std::map<std::string, DispersionRelation<double>> &
                dispersion_relations) {
            std::vector<SphericalECGDispersionContext> contexts;
            for (const auto &pair : dispersion_relations) {
                contexts.emplace_back(space_dimension, particles, pair.second);
            }
            return contexts;
        }

        static std::vector<SphericalECGConfiningContext>
        extract_confining_contexts() {
            std::vector<SphericalECGConfiningContext> contexts;
            // TODO
            return contexts;
        }

        static std::vector<SphericalECGPairwiseContext>
        extract_pairwise_contexts(
                long long int space_dimension,
                const std::vector<Particle<double>> &particles,
                const std::map<std::string, PairwisePotential<double>> &
                pairwise_potentials) {
            std::vector<SphericalECGPairwiseContext> contexts;
            for (const auto &pair : pairwise_potentials) {
                contexts.emplace_back(space_dimension, particles, pair.second);
            }
            return contexts;
        }

        static std::vector<bool> extract_permutation_signatures(
                const std::vector<Particle<double>> &particles,
                const std::vector<std::vector<std::size_t>> &permutations) {
            std::vector<bool> signatures;
            for (const auto &permutation : permutations) {
                const std::size_t signature =
                        dznl::count_changes(particles, permutation) / 2 +
                        dznl::count_inversions(permutation);
                signatures.push_back(static_cast<bool>(signature % 2));
            }
            return signatures;
        }

        const long long int space_dimension;
        const std::size_t num_particles;
        const std::size_t num_parameters;
        packed_determinant_inverse_function<double> packed_determinant_inverse;
        std::vector<double> a;
        std::vector<double> b;
        std::vector<double> c;
        std::vector<double> d;
        const std::vector<std::vector<std::size_t>> allowed_permutations;
        const std::size_t num_permutations;
        const std::vector<bool> permutation_signatures;
        const std::vector<SphericalECGDispersionContext> dispersion_contexts;
        const std::vector<SphericalECGConfiningContext> confining_contexts;
        const std::vector<SphericalECGPairwiseContext> pairwise_contexts;

    public:

        explicit SphericalECGOverlapContext(
                long long int space_dimension,
                const std::vector<Particle<double>> &particles,
                const std::map<std::string, DispersionRelation<double>> &
                dispersion_relations,
                const std::map<std::string, ConfiningPotential<double>> &,
                const std::map<std::string, PairwisePotential<double>> &
                pairwise_potentials)
                : space_dimension(space_dimension),
                  num_particles(particles.size()),
                  num_parameters(particles.size() * (particles.size() + 1) / 2),
                  packed_determinant_inverse(
                          PACKED_DETERMINANT_INVERSE<double>[num_particles -
                                                             1]),
                  a(num_parameters),
                  b(num_parameters),
                  c(num_parameters),
                  d(num_parameters),
                  allowed_permutations(dznl::invariant_permutations(particles)),
                  num_permutations(allowed_permutations.size()),
                  permutation_signatures(extract_permutation_signatures(
                          particles, allowed_permutations)),
                  dispersion_contexts(extract_dispersion_contexts(
                          space_dimension, particles, dispersion_relations)),
                  confining_contexts(extract_confining_contexts()),
                  pairwise_contexts(extract_pairwise_contexts(
                          space_dimension, particles, pairwise_potentials)) {}

        void evaluate_kernel(
                double &__restrict__ overlap_kernel,
                double &__restrict__ hamiltonian_kernel) {
            for (std::size_t k = 0; k < num_parameters; ++k) {
                c[k] = a[k] + b[k];
            }
            overlap_kernel = half_inverse_pow(
                    packed_determinant_inverse(c.data(), d.data()),
                    space_dimension);
            hamiltonian_kernel = 0.0;
            for (const auto &dispersion_context : dispersion_contexts) {
                hamiltonian_kernel += dispersion_context.evaluate_kernel(
                        a.data(), b.data(), d.data());
            }
            for (const auto &confining_context : confining_contexts) {
                // TODO
                hamiltonian_kernel += confining_context.evaluate_kernel();
            }
            for (const auto &pairwise_context : pairwise_contexts) {
                hamiltonian_kernel += pairwise_context.evaluate_kernel(
                        d.data());
            }
            hamiltonian_kernel *= overlap_kernel;
        }

        static std::size_t triangular_index(
                std::size_t i, std::size_t j) noexcept {
            const std::size_t lo = std::min(i, j);
            const std::size_t hi = std::max(i, j);
            return hi * (hi + 1) / 2 + lo;
        }

        void evaluate_matrix_elements(
                double &__restrict__ overlap_element,
                double &__restrict__ hamiltonian_element,
                const double *aa, const double *bb) {
            overlap_element = hamiltonian_element = 0.0;
            for (std::size_t p = 0; p < num_permutations; ++p) {
                const auto &perm_a = allowed_permutations[p];
                for (std::size_t q = 0; q < num_permutations; ++q) {
                    const auto &perm_b = allowed_permutations[q];
                    for (std::size_t i = 0, k = 0; i < num_particles; ++i) {
                        for (std::size_t j = 0; j <= i; ++j, ++k) {
                            a[k] = aa[triangular_index(perm_a[i], perm_a[j])];
                            b[k] = bb[triangular_index(perm_b[i], perm_b[j])];
                        }
                    }
                    double overlap_kernel, hamiltonian_kernel;
                    evaluate_kernel(overlap_kernel, hamiltonian_kernel);
                    if (permutation_signatures[p] ^ permutation_signatures[q]) {
                        overlap_element -= overlap_kernel;
                        hamiltonian_element -= hamiltonian_kernel;
                    } else {
                        overlap_element += overlap_kernel;
                        hamiltonian_element += hamiltonian_kernel;
                    }
                }
            }
        }

    }; // class SphericalECGOverlapContext

} // namespace zsvm

#endif // ZSVM_SPHERICAL_ECG_OVERLAP_CONTEXT_HPP_INCLUDED
