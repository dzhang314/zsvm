#ifndef ZSVM_REAL_SOLVER_HPP_INCLUDED
#define ZSVM_REAL_SOLVER_HPP_INCLUDED

// C++ standard library headers
#include <cmath> // for std::tgamma
#include <cstddef> // for std::size_t
#include <iostream> // for std::cout
#include <stdexcept> // for std::invalid_argument
#include <vector> // for std::vector

// Eigen linear algebra library headers
#include <Eigen/Core>
#include <Eigen/LU>

// Project-specific headers
#include "Particle.hpp"
#include "Permutation.hpp"
#include "JacobiCoordinates.hpp"

namespace zsvm {

    class RealSolver {

    private: // =============================================== MEMBER VARIABLES

        const int space_dimension;
        const std::size_t num_particles;
        const std::size_t num_pairs;

        std::vector<int> particle_types;
        std::vector<double> masses;
        std::vector<double> charges;
        std::vector<Spin> spins;

        Eigen::MatrixXd inverse_mass_matrix;
        Eigen::MatrixXd pairwise_weight_vectors;
        std::vector<Eigen::MatrixXd> pairwise_weight_matrices;
        std::vector<double> pairwise_charge_products;
        std::vector<std::vector<std::size_t>> allowed_permutations;
        Eigen::MatrixXd permutation_sign_matrix;
        std::vector<Eigen::MatrixXd> permutation_matrices;

    public: // ===================================================== CONSTRUCTOR

        explicit RealSolver(const std::vector<Particle> &particles,
                            int space_dimension)
                : space_dimension(space_dimension),
                  num_particles(particles.size()),
                  num_pairs(particles.size() * (particles.size() - 1) / 2),
                  particle_types(particles.size()),
                  masses(particles.size()),
                  charges(particles.size()),
                  spins(particles.size()),
                  pairwise_charge_products(
                          particles.size() * (particles.size() - 1) / 2) {
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
            pairwise_weight_vectors = jaco::pairwise_weights(masses);
            for (std::size_t k = 0; k < num_pairs; ++k) {
                pairwise_weight_matrices.push_back(
                        pairwise_weight_vectors.col(k) *
                        pairwise_weight_vectors.col(k).transpose());
            }
            for (std::size_t i = 0, k = 0; i < num_particles - 1; ++i) {
                for (std::size_t j = i + 1; j < num_particles; ++j, ++k) {
                    pairwise_charge_products[k] = charges[i] * charges[j];
                }
            }
            allowed_permutations = dznl::invariant_permutations(particle_types);
            permutation_sign_matrix.setConstant(
                    allowed_permutations.size(), allowed_permutations.size(),
                    0.0);
            for (std::size_t i = 0; i < allowed_permutations.size(); ++i) {
                for (std::size_t j = 0; j < allowed_permutations.size(); ++j) {
                    const std::size_t signature =
                            dznl::count_changes(
                                    spins, allowed_permutations[i]) / 2 +
                            dznl::count_changes(
                                    spins, allowed_permutations[j]) / 2 +
                            dznl::count_inversions(allowed_permutations[i]) +
                            dznl::count_inversions(allowed_permutations[j]);
                    permutation_sign_matrix(i, j) = (signature % 2 == 0)
                                                    ? +1.0 : -1.0;
                }
            }
            permutation_matrices =
                    jaco::permutation_matrices(masses, allowed_permutations);
        }

    public: // =================================== MATRIX ELEMENT HELPER METHODS

        Eigen::MatrixXd gaussian_parameter_matrix(
                const std::vector<double> &correlation_coefficients) const {
            if (correlation_coefficients.size() != num_pairs) {
                throw std::invalid_argument(
                        "RealSolver::gaussian_parameter_matrix received "
                                "vector with incorrect number of "
                                "correlation coefficients");
            }
            Eigen::MatrixXd a =
                    correlation_coefficients[0] * pairwise_weight_matrices[0];
            for (std::size_t k = 1; k < num_pairs; ++k) {
                a += correlation_coefficients[k] * pairwise_weight_matrices[k];
            }
            return a;
        }

        double overlap_kernel(const Eigen::MatrixXd &a,
                              const Eigen::MatrixXd &b) const {
            return std::sqrt(std::pow((a + b).determinant(), -space_dimension));
        }

        double kinetic_kernel(const Eigen::MatrixXd &a,
                              const Eigen::MatrixXd &b) const {
            return (0.5 * space_dimension) * overlap_kernel(a, b) *
                   (inverse_mass_matrix * a * (a + b).inverse() * b).trace();
        }

        double coulomb_kernel(const Eigen::MatrixXd &a,
                              const Eigen::MatrixXd &b) const {
            const Eigen::MatrixXd inv = (a + b).inverse();
            const double dimension_factor =
                    std::tgamma(0.5 * (space_dimension - 1.0)) /
                    std::tgamma(0.5 * space_dimension);
            double result = 0.0;
            for (std::size_t k = 0; k < num_pairs; ++k) {
                result += pairwise_charge_products[k] / std::sqrt(
                        2.0 * pairwise_weight_vectors.col(k).dot(
                                inv * pairwise_weight_vectors.col(k)));
            }
            return overlap_kernel(a, b) * dimension_factor * result;
        }

    public: // ========================================== MATRIX ELEMENT METHODS

        double overlap_matrix_element(const Eigen::MatrixXd &a,
                                      const Eigen::MatrixXd &b) const {
            double result = 0.0;
            for (std::size_t i = 0; i < allowed_permutations.size(); ++i) {
                for (std::size_t j = 0; j < allowed_permutations.size(); ++j) {
                    result += permutation_sign_matrix(i, j) * overlap_kernel(
                            permutation_matrices[i].transpose() *
                            a * permutation_matrices[i],
                            permutation_matrices[j].transpose() *
                            b * permutation_matrices[j]);
                }
            }
            return result;
        }

        double kinetic_matrix_element(const Eigen::MatrixXd &a,
                                      const Eigen::MatrixXd &b) const {
            double result = 0.0;
            for (std::size_t i = 0; i < allowed_permutations.size(); ++i) {
                for (std::size_t j = 0; j < allowed_permutations.size(); ++j) {
                    result += permutation_sign_matrix(i, j) * kinetic_kernel(
                            permutation_matrices[i].transpose() *
                            a * permutation_matrices[i],
                            permutation_matrices[j].transpose() *
                            b * permutation_matrices[j]);
                }
            }
            return result;
        }

        double coulomb_matrix_element(const Eigen::MatrixXd &a,
                                      const Eigen::MatrixXd &b) const {
            double result = 0.0;
            for (std::size_t i = 0; i < allowed_permutations.size(); ++i) {
                for (std::size_t j = 0; j < allowed_permutations.size(); ++j) {
                    result += permutation_sign_matrix(i, j) * coulomb_kernel(
                            permutation_matrices[i].transpose() *
                            a * permutation_matrices[i],
                            permutation_matrices[j].transpose() *
                            b * permutation_matrices[j]);
                }
            }
            return result;
        }

    }; // class RealSolver

} // namespace zsvm

#endif // ZSVM_REAL_SOLVER_HPP_INCLUDED
