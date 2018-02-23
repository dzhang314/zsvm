#ifndef ZSVM_JACOBI_COORDINATES_HPP_INCLUDED
#define ZSVM_JACOBI_COORDINATES_HPP_INCLUDED

// C++ standard library headers
#include <vector>

// Eigen linear algebra library headers
#include <Eigen/Core>

namespace jaco {

    Eigen::MatrixXd transformation_matrix(
            const std::vector<double> &masses) {
        const std::size_t n = masses.size();
        std::vector<double> mass_sums(n);
        if (n > 0) { mass_sums[0] = masses[0]; }
        for (std::size_t i = 1; i < n; ++i) {
            mass_sums[i] = mass_sums[i - 1] + masses[i];
        }
        Eigen::MatrixXd u(n, n);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                if (i >= j) {
                    u(i, j) = masses[j] / mass_sums[i];
                } else if (i == j - 1) {
                    u(i, j) = -1.0;
                } else {
                    u(i, j) = 0.0;
                }
            }
        }
        return u;
    }

    Eigen::MatrixXd transformation_matrix_inverse(
            const std::vector<double> &masses) {
        const std::size_t n = masses.size();
        std::vector<double> mass_sums(n);
        if (n > 0) { mass_sums[0] = masses[0]; }
        for (std::size_t i = 1; i < n; ++i) {
            mass_sums[i] = mass_sums[i - 1] + masses[i];
        }
        Eigen::MatrixXd v(n, n);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                if (j + 1 == n) {
                    v(i, j) = 1.0;
                } else if (i <= j) {
                    v(i, j) = masses[j + 1] / mass_sums[j + 1];
                } else if (i == j + 1) {
                    v(i, j) = -mass_sums[j] / mass_sums[i];
                } else {
                    v(i, j) = 0.0;
                }
            }
        }
        return v;
    }

} // namespace jaco

#endif // ZSVM_JACOBI_COORDINATES_HPP_INCLUDED
