#ifndef ZSVM_JACOBI_COORDINATES_HPP_INCLUDED
#define ZSVM_JACOBI_COORDINATES_HPP_INCLUDED

// C++ standard library headers
#include <cstddef> // for std::size_t
#include <vector> // for std::vector

// Eigen linear algebra library headers
#include <Eigen/Core> // for Eigen::MatrixXd

namespace jaco {

    Eigen::MatrixXd transformation_matrix(
            const std::vector<double> &masses);

    Eigen::MatrixXd transformation_matrix_inverse(
            const std::vector<double> &masses);

    Eigen::MatrixXd inverse_mass_matrix(
            const std::vector<double> &masses);

    Eigen::MatrixXd reduced_inverse_mass_matrix(
            const std::vector<double> &masses);

    Eigen::MatrixXd pairwise_weights(
            const std::vector<double> &masses);

    std::vector<Eigen::MatrixXd> permutation_matrices(
            const std::vector<double> &masses,
            const std::vector<std::vector<std::size_t>> &permutations);

} // namespace jaco

#endif // ZSVM_JACOBI_COORDINATES_HPP_INCLUDED
