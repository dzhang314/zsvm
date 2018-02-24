#ifndef ZSVM_JACOBI_COORDINATES_HPP_INCLUDED
#define ZSVM_JACOBI_COORDINATES_HPP_INCLUDED

// C++ standard library headers
#include <stdexcept>
#include <vector>

// Eigen linear algebra library headers
#include <Eigen/Core>

namespace jaco {

    Eigen::MatrixXd transformation_matrix(const std::vector<double> &masses) {
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

    Eigen::MatrixXd inverse_mass_matrix(const std::vector<double> &masses) {
        const std::size_t n = masses.size();
        const Eigen::MatrixXd u = transformation_matrix(masses);
        Eigen::MatrixXd inverse_masses(n, n);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                if (i == j) {
                    inverse_masses(i, j) = 1.0 / masses[i];
                } else {
                    inverse_masses(i, j) = 0.0;
                }
            }
        }
        return u * inverse_masses * u.transpose();
    }

    Eigen::MatrixXd reduced_inverse_mass_matrix(
            const std::vector<double> &masses) {
        const std::size_t n = masses.size();
        if (n == 0) {
            throw std::invalid_argument(
                    "reduced_inverse_mass_matrix "
                            "received empty vector of masses");
        }
        const Eigen::MatrixXd m = inverse_mass_matrix(masses);
        Eigen::MatrixXd r(n - 1, n - 1);
        for (std::size_t i = 0; i < n - 1; ++i) {
            for (std::size_t j = 0; j < n - 1; ++j) {
                r(i, j) = m(i, j);
            }
        }
        return r;
    }

    Eigen::MatrixXd pairwise_weights(const std::vector<double> &masses) {
        const std::size_t n = masses.size();
        if (n == 0) {
            throw std::invalid_argument(
                    "pairwise_weights received empty vector of masses");
        }
        const std::size_t num_pairs = n * (n - 1) / 2;
        const Eigen::MatrixXd v = transformation_matrix_inverse(masses);
        Eigen::MatrixXd result(n - 1, num_pairs);
        for (std::size_t i = 0, k = 0; i < n - 1; ++i) {
            for (std::size_t j = i + 1; j < n; ++j, ++k) {
                for (std::size_t m = 0; m < n - 1; ++m) {
                    result(m, k) = v(i, m) - v(j, m);
                }
            }
        }
        return result;
    }

    std::vector<Eigen::MatrixXd> permutation_matrices(
            const std::vector<double> &masses,
            const std::vector<std::vector<std::size_t>> &permutations) {
        const std::size_t n = masses.size();
        if (n == 0) {
            throw std::invalid_argument(
                    "permutation_matrices "
                            "received empty vector of masses");
        }
        const Eigen::MatrixXd u = transformation_matrix(masses);
        const Eigen::MatrixXd v = transformation_matrix_inverse(masses);
        std::vector<Eigen::MatrixXd> result;
        for (const auto &permutation : permutations) {
            Eigen::MatrixXd w(n, n);
            for (std::size_t i = 0; i < n; ++i) {
                for (std::size_t j = 0; j < n; ++j) {
                    w(i, j) = v(permutation[i], j);
                }
            }
            const Eigen::MatrixXd x = u * w;
            result.emplace_back(n - 1, n - 1);
            for (std::size_t i = 0; i < n - 1; ++i) {
                for (std::size_t j = 0; j < n - 1; ++j) {
                    result.back()(i, j) = x(i, j);
                }
            }
        }
        return result;
    }

} // namespace jaco

#endif // ZSVM_JACOBI_COORDINATES_HPP_INCLUDED
