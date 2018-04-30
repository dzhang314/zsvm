#ifndef ZSVM_SPHERICAL_ECG_JACOBI_CONTEXT_HPP
#define ZSVM_SPHERICAL_ECG_JACOBI_CONTEXT_HPP

// C++ standard library headers
#include <cstddef> // for std::size_t
#include <vector>

// Eigen linear algebra library headers
#include <Eigen/Core>

// Project-specific headers
#include "JacobiCoordinates.hpp"
#include "PackedLinearAlgebra.hpp"
#include "Particle.hpp"
#include "Permutation.hpp"

namespace zsvm {

    template <typename T>
    class SphericalECGJacobiContext {

    private: // =============================================== MEMBER VARIABLES

        typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> MatrixXT;

        const std::size_t num_particles;
        const std::size_t num_pairs;
        const std::size_t num_permutations;
        const std::size_t matrix_size;
        const long long int space_dimension;
        const T dimension_factor;
        const T kinetic_factor;

        const std::vector<T> inverse_masses;
        const std::vector<T> weight_vectors;
        const std::vector<T> weight_matrices;
        const std::vector<T> charge_products;
        const std::vector<T> permutation_signs;
        const std::vector<T> permutation_matrices;

        std::vector<T> vax;
        std::vector<T> vbx;
        std::vector<T> vcx;
        std::vector<T> vdx;

        const packed_determinant_inverse_function<T> packed_determinant_inverse;
        const packed_kinetic_trace_function<T> packed_kinetic_trace;
        const packed_quadratic_form_function<T> packed_quadratic_form;
        const packed_permutation_conjugate_function<T> packed_permutation_conjugate;

    private: // ============================================ FACTORY CONSTRUCTOR

        T gamma(const T &x) {
            using std::tgamma;
            return tgamma(x);
        }

        explicit SphericalECGJacobiContext(
                const std::vector<T> &charges,
                std::size_t num_permutations,
                long long int space_dimension,
                const MatrixXT &inverse_mass_matrix,
                const MatrixXT &pairwise_weight_vectors,
                const std::vector<T> &permutation_sign_vector,
                const std::vector<MatrixXT> &permutation_matrix_vector)
                : num_particles(charges.size()),
                  num_pairs(num_particles * (num_particles - 1) / 2),
                  num_permutations(num_permutations),
                  matrix_size((num_particles - 1) * (num_particles - 1)),
                  space_dimension(space_dimension),
                  dimension_factor(gamma((space_dimension - static_cast<T>(1)) /
                                         static_cast<T>(2)) /
                                   gamma(space_dimension / static_cast<T>(2))),
                  kinetic_factor(space_dimension / (2 * dimension_factor)),
                  inverse_masses(num_particles - 1),
                  weight_vectors(num_pairs * (num_particles - 1)),
                  weight_matrices(num_pairs * num_pairs),
                  charge_products(num_pairs),
                  permutation_signs(num_permutations),
                  permutation_matrices(num_permutations * matrix_size),
                  vax(num_pairs),
                  vbx(num_pairs),
                  vcx(num_pairs),
                  vdx(num_pairs),
                  packed_determinant_inverse(
                          PackedDeterminantInverse<T>::FUNC[
                                  num_particles - 2]),
                  packed_kinetic_trace(
                          PackedKineticTrace<T>::FUNC[
                                  num_particles - 2]),
                  packed_quadratic_form(
                          PackedQuadraticForm<T>::FUNC[
                                  num_particles - 2]),
                  packed_permutation_conjugate(
                          PackedPermutationConjugate<T>::FUNC[
                                  num_particles - 2]) {
            auto inverse_mass_pointer = const_cast<T *>(inverse_masses.data());
            for (std::size_t i = 0; i < num_particles - 1; ++i) {
                inverse_mass_pointer[i] = inverse_mass_matrix(i, i);
            }
            auto weight_vector_pointer = const_cast<T *>(weight_vectors.data());
            for (std::size_t p = 0, k = 0; p < num_pairs; ++p) {
                for (std::size_t i = 0; i < num_particles - 1; ++i, ++k) {
                    weight_vector_pointer[k] = pairwise_weight_vectors(i, p);
                }
            }
            auto weight_matrix_pointer = const_cast<T *>(weight_matrices.data());
            for (std::size_t p = 0, k = 0; p < num_pairs; ++p) {
                for (std::size_t i = 0; i < num_particles - 1; ++i) {
                    for (std::size_t j = 0; j <= i; ++j, ++k) {
                        weight_matrix_pointer[k] =
                                pairwise_weight_vectors(i, p) *
                                pairwise_weight_vectors(j, p);
                    }
                }
            }
            auto charge_product_pointer = const_cast<T *>(charge_products.data());
            for (std::size_t i = 0, k = 0; i < num_particles - 1; ++i) {
                for (std::size_t j = i + 1; j < num_particles; ++j, ++k) {
                    charge_product_pointer[k] = charges[i] * charges[j];
                }
            }
            auto permutation_sign_pointer = const_cast<T *>(permutation_signs.data());
            for (std::size_t p = 0; p < num_permutations; ++p) {
                permutation_sign_pointer[p] = permutation_sign_vector[p];
            }
            auto permutation_matrix_pointer = const_cast<T *>(permutation_matrices.data());
            for (std::size_t p = 0, k = 0; p < num_permutations; ++p) {
                for (std::size_t i = 0; i < num_particles - 1; ++i) {
                    for (std::size_t j = 0; j < num_particles - 1; ++j, ++k) {
                        permutation_matrix_pointer[k] =
                                permutation_matrix_vector[p](j, i);
                    }
                }
            }
        }

    public: // =========================================== STATIC FACTORY METHOD

        static SphericalECGJacobiContext create(
                const std::vector<zsvm::Particle<T>> &particles,
                const std::string &mass_carrier,
                const std::string &charge_carrier,
                long long int space_dimension) {
            const std::size_t num_particles = particles.size();
            if (num_particles < 2) {
                throw std::invalid_argument(
                        "Attempted to construct SphericalECGJacobiContext "
                        "with fewer than 2 particles");
            }
            std::vector<T> masses(num_particles);
            std::vector<T> charges(num_particles);
            for (std::size_t i = 0; i < num_particles; ++i) {
                masses[i] = particles[i].carriers.at(mass_carrier);
                charges[i] = particles[i].carriers.at(charge_carrier);
            }
            const std::vector<std::vector<std::size_t>> allowed_permutations =
                    dznl::invariant_permutations(particles);
            std::vector<T> permutation_signs;
            for (const auto &permutation : allowed_permutations) {
                const std::size_t signature =
                        dznl::count_changes(particles, permutation) / 2 +
                        dznl::count_inversions(permutation);
                permutation_signs.push_back((signature % 2 == 0) ? +1 : -1);
            }
            return SphericalECGJacobiContext(
                    charges,
                    allowed_permutations.size(),
                    space_dimension,
                    jaco::reduced_inverse_mass_matrix(masses),
                    jaco::pairwise_weights(masses),
                    permutation_signs,
                    jaco::permutation_matrices(masses, allowed_permutations));
        }

    private: // ================================== MATRIX ELEMENT HELPER METHODS

        void matrix_element_kernel(
                T &__restrict__ overlap_kernel,
                T &__restrict__ hamiltonian_kernel) {
            using std::sqrt;
            T *__restrict__ const ax = vax.data();
            T *__restrict__ const bx = vbx.data();
            T *__restrict__ const cx = vcx.data();
            T *__restrict__ const dx = vdx.data();
            for (std::size_t i = 0; i < num_pairs; ++i) {
                cx[i] = ax[i] + bx[i];
            }
            overlap_kernel = half_inverse_pow(
                    packed_determinant_inverse(cx, dx), space_dimension);
            hamiltonian_kernel = kinetic_factor * packed_kinetic_trace(
                    ax, bx, dx, inverse_masses.data());
            for (std::size_t k = 0; k < num_pairs; ++k) {
                const T alpha = packed_quadratic_form(
                        dx, weight_vectors.data() + (num_particles - 1) * k);
                hamiltonian_kernel += charge_products[k] / sqrt(2 * alpha);
            }
            hamiltonian_kernel *= dimension_factor * overlap_kernel;
        }

    public: // ========================================== MATRIX ELEMENT METHODS

        void gaussian_parameter_matrix(
                const T *__restrict__ correlation_coefficients,
                T *__restrict__ result) const {
            using std::exp;
            for (std::size_t p = 0; p < num_pairs; ++p) { result[p] = 0; }
            for (std::size_t p = 0, k = 0; p < num_pairs; ++p) {
                const T c = exp(correlation_coefficients[p]);
                for (std::size_t q = 0; q < num_pairs; ++q, ++k) {
                    result[q] += c * weight_matrices[k];
                }
            }
        }

        void evaluate_matrix_elements(
                T &__restrict__ overlap_element,
                T &__restrict__ hamiltonian_element,
                const T *a, const T *b) {
            T *__restrict__ const ax = vax.data();
            T *__restrict__ const bx = vbx.data();
            overlap_element = hamiltonian_element = 0;
            T overlap_kernel, hamiltonian_kernel;
            for (std::size_t i = 0; i < num_permutations; ++i) {
                for (std::size_t j = 0; j < num_permutations; ++j) {
                    const T sign = permutation_signs[i] * permutation_signs[j];
                    packed_permutation_conjugate(
                            a, permutation_matrices.data() + i * matrix_size,
                            ax);
                    packed_permutation_conjugate(
                            b, permutation_matrices.data() + j * matrix_size,
                            bx);
                    matrix_element_kernel(overlap_kernel, hamiltonian_kernel);
                    overlap_element += sign * overlap_kernel;
                    hamiltonian_element += sign * hamiltonian_kernel;
                }
            }
        }

    }; // class SphericalECGJacobiContext

} // namespace zsvm

#endif // ZSVM_SPHERICAL_ECG_JACOBI_CONTEXT_HPP
