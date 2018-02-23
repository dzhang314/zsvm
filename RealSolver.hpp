#ifndef ZSVM_REAL_SOLVER_HPP_INCLUDED
#define ZSVM_REAL_SOLVER_HPP_INCLUDED

// C++ standard library headers
#include <cstdint>
#include <iostream>
#include <vector>

// Eigen linear algebra library headers
#include <Eigen/Core>

// Project-specific headers
#include "Particle.hpp"
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
                                "with fewer than 2 particles.");
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
            const Eigen::MatrixXd u =
                    jaco::transformation_matrix(masses);
            const Eigen::MatrixXd v =
                    jaco::transformation_matrix_inverse(masses);

        }

    }; // class RealSolver

} // namespace zsvm

#endif // ZSVM_REAL_SOLVER_HPP_INCLUDED
