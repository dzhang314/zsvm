#ifndef ZSVM_JACOBI_COORDINATES_HPP_INCLUDED
#define ZSVM_JACOBI_COORDINATES_HPP_INCLUDED

// C++ standard library headers
#include <cstddef> // for std::size_t
#include <vector>

// Eigen linear algebra library headers
#include <Eigen/Core>

namespace jaco {

    template <typename T>
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
    transformation_matrix(const std::vector<T> &masses) {
        const std::size_t n = masses.size();
        std::vector<T> mass_sums(n);
        if (n > 0) { mass_sums[0] = masses[0]; }
        for (std::size_t i = 1; i < n; ++i) {
            mass_sums[i] = mass_sums[i - 1] + masses[i];
        }
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> u(n, n);
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

    template <typename T>
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
    transformation_matrix_inverse(const std::vector<T> &masses) {
        const std::size_t n = masses.size();
        std::vector<T> mass_sums(n);
        if (n > 0) { mass_sums[0] = masses[0]; }
        for (std::size_t i = 1; i < n; ++i) {
            mass_sums[i] = mass_sums[i - 1] + masses[i];
        }
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> v(n, n);
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

    template <typename T>
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
    inverse_mass_matrix(const std::vector<T> &masses) {
        const std::size_t n = masses.size();
        const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> u =
                transformation_matrix(masses);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> inverse_masses(n, n);
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

    template <typename T>
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
    reduced_inverse_mass_matrix(const std::vector<T> &masses) {
        const std::size_t n = masses.size();
        if (n == 0) {
            throw std::invalid_argument(
                    "reduced_inverse_mass_matrix "
                    "received empty vector of masses");
        }
        const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> m =
                inverse_mass_matrix(masses);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> r(n - 1, n - 1);
        for (std::size_t i = 0; i < n - 1; ++i) {
            for (std::size_t j = 0; j < n - 1; ++j) {
                r(i, j) = m(i, j);
            }
        }
        return r;
    }

    template <typename T>
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
    pairwise_weights(const std::vector<T> &masses) {
        const std::size_t n = masses.size();
        if (n == 0) {
            throw std::invalid_argument(
                    "pairwise_weights received empty vector of masses");
        }
        const std::size_t num_pairs = n * (n - 1) / 2;
        const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> v =
                transformation_matrix_inverse(masses);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> w(n - 1, num_pairs);
        for (std::size_t i = 0, k = 0; i < n - 1; ++i) {
            for (std::size_t j = i + 1; j < n; ++j, ++k) {
                for (std::size_t m = 0; m < n - 1; ++m) {
                    w(m, k) = v(i, m) - v(j, m);
                }
            }
        }
        return w;
    }

    template <typename T>
    std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>
    permutation_matrices(
            const std::vector<T> &masses,
            const std::vector<std::vector<std::size_t>> &permutations) {
        const std::size_t n = masses.size();
        if (n == 0) {
            throw std::invalid_argument(
                    "permutation_matrices received empty vector of masses");
        }
        const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> u =
                transformation_matrix(masses);
        const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> v =
                transformation_matrix_inverse(masses);
        std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> result;
        for (const auto &permutation : permutations) {
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> w(n, n);
            for (std::size_t i = 0; i < n; ++i) {
                for (std::size_t j = 0; j < n; ++j) {
                    w(i, j) = v(permutation[i], j);
                }
            }
            const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> x = u * w;
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
