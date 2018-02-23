#ifndef ZSVM_REAL_SOLVER_HPP_INCLUDED
#define ZSVM_REAL_SOLVER_HPP_INCLUDED

// C++ standard library headers
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <vector>

// Eigen linear algebra library headers
#include <Eigen/Core>

// Project-specific headers
#include "Particle.hpp"
#include "Permutation.hpp"
#include "JacobiCoordinates.hpp"

namespace zsvm {

    class RealSolver {

    private: // =============================================== MEMBER VARIABLES

        const std::size_t num_particles;
        const std::size_t num_pairs;

        std::vector<int> particle_types;
        std::vector<double> masses;
        std::vector<double> charges;
        std::vector<Spin> spins;

        Eigen::MatrixXd inverse_mass_matrix;
        Eigen::MatrixXd pairwise_weights;
        std::vector<std::vector<std::size_t>> allowed_permutations;
        std::vector<Eigen::MatrixXd> permutation_matrices;

    public: // ===================================================== CONSTRUCTOR

        explicit RealSolver(const std::vector<Particle> &particles)
                : num_particles(particles.size()),
                  num_pairs(particles.size() * (particles.size() - 1) / 2),
                  particle_types(particles.size()),
                  masses(particles.size()),
                  charges(particles.size()),
                  spins(particles.size()) {
            if (num_particles < 2) {
                throw std::invalid_argument(
                        "Attempted to construct SVMSolver "
                                "with fewer than 2 particles");
            }
            std::cout << "Constructing SVMSolver with "
                      << num_particles << " particles.\n";
            for (std::size_t i = 0; i < num_particles; ++i) {
                std::cout << "Particle " << i << ": " << particles[i] << "\n";
                particle_types[i] = particles[i].type;
                masses[i] = particles[i].mass;
                charges[i] = particles[i].charge;
                spins[i] = particles[i].spin;
            }
            inverse_mass_matrix = jaco::reduced_inverse_mass_matrix(masses);
            std::cout << inverse_mass_matrix << std::endl;
            pairwise_weights = jaco::pairwise_weights(masses);
            std::cout << pairwise_weights << std::endl;
            allowed_permutations = dznl::invariant_permutations(particle_types);
            std::cout << "Order of particle exchange symmetry group: "
                      << allowed_permutations.size() << '\n';
            for (const auto &permutation : allowed_permutations) {
                for (const auto &index : permutation) {
                    std::cout << index;
                }
                std::cout << '\n';
            }
            permutation_matrices =
                    jaco::permutation_matrices(masses, allowed_permutations);
            for (const auto &matrix : permutation_matrices) {
                std::cout << matrix << '\n';
            }
        }

    }; // class RealSolver

} // namespace zsvm

#endif // ZSVM_REAL_SOLVER_HPP_INCLUDED
