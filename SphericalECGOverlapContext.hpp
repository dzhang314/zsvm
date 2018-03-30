#ifndef ZSVM_SPHERICAL_ECG_OVERLAP_CONTEXT_HPP_INCLUDED
#define ZSVM_SPHERICAL_ECG_OVERLAP_CONTEXT_HPP_INCLUDED

#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include "Particle.hpp"
#include "Permutation.hpp"
#include "SphericalECGDispersionContext.hpp"
#include "SphericalECGConfiningContext.hpp"
#include "SphericalECGPairwiseContext.hpp"

namespace zsvm {

    class SphericalECGOverlapContext {

        static std::vector<SphericalECGDispersionContext>
        extract_dispersion_contexts(
                const std::vector<Particle> &particles,
                const std::map<std::string, DispersionRelation> &dispersion_relations) {
            std::vector<SphericalECGDispersionContext> result;
            for (const auto &pair : dispersion_relations) {
                result.emplace_back(particles, pair.second);
            }
            return result;
        }

        static std::vector<SphericalECGConfiningContext>
        extract_confining_contexts() {
            std::vector<SphericalECGConfiningContext> result;
            // TODO
            return result;
        }

        static std::vector<SphericalECGPairwiseContext>
        extract_pairwise_contexts(
                long long int space_dimension,
                const std::vector<Particle> &particles,
                const std::map<std::string, PairwisePotential> &pairwise_potentials) {
            std::vector<SphericalECGPairwiseContext> result;
            for (const auto &pair : pairwise_potentials) {
                result.emplace_back(space_dimension, particles, pair.second);
            }
            return result;
        }

        const std::size_t num_particles;
        const std::size_t num_parameters;
        const std::vector<std::vector<std::size_t>> allowed_permutations;
        const std::vector<SphericalECGDispersionContext> dispersion_contexts;
        const std::vector<SphericalECGConfiningContext> confining_contexts;
        const std::vector<SphericalECGPairwiseContext> pairwise_contexts;

        explicit SphericalECGOverlapContext(
                long long int space_dimension,
                const std::vector<Particle> &particles,
                const std::map<std::string, DispersionRelation> &dispersion_relations,
                const std::map<std::string, ConfiningPotential> &confining_potentials,
                const std::map<std::string, PairwisePotential> &pairwise_potentials)
                : num_particles(particles.size()),
                  num_parameters(particles.size() * (particles.size() + 1) / 2),
                  allowed_permutations(dznl::invariant_permutations(particles)),
                  dispersion_contexts(extract_dispersion_contexts(
                          particles, dispersion_relations)),
                  confining_contexts(extract_confining_contexts()),
                  pairwise_contexts(extract_pairwise_contexts(
                          space_dimension, particles, pairwise_potentials)) {}

        void evaluate_kernel(
                double &__restrict__ overlap_kernel,
                double &__restrict__ hamiltonian_kernel) {
            overlap_kernel = 0.0;
            hamiltonian_kernel = 0.0;
        }

        void evaluate_matrix_elements(
                double &__restrict__ overlap_element,
                double &__restrict__ hamiltonian_element) {
            overlap_element = 0.0;
            hamiltonian_element = 0.0;
        }

    }; // class SphericalECGOverlapContext

} // namespace zsvm

#endif // ZSVM_SPHERICAL_ECG_OVERLAP_CONTEXT_HPP_INCLUDED
