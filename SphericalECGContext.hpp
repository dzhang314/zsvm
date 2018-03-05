#ifndef ZSVM_SPHERICAL_ECG_CONTEXT_HPP_INCLUDED
#define ZSVM_SPHERICAL_ECG_CONTEXT_HPP_INCLUDED

// C++ standard library headers
#include <cstddef> // for std::size_t
#include <random> // for std::mt19937_64
#include <vector> // for std::vector

// Eigen linear algebra library headers
#include <Eigen/Core> // for Eigen::MatrixXd

// Project-specific headers
#include "Particle.hpp" // for zsvm::Particle

#define ZSVM_SPHERICAL_ECG_CONTEXT_TIMING_ENABLED 1

namespace zsvm {

    class SphericalECGContext {

    private: // =============================================== MEMBER VARIABLES

        const std::size_t num_particles;
        const std::size_t num_pairs;
        const std::size_t num_permutations;
        const std::size_t matrix_size;
        const int space_dimension;
        const double dimension_factor;
        const double kinetic_factor;

        const double *__restrict__ const inverse_masses;
        const double *__restrict__ const weight_vectors;
        const double *__restrict__ const weight_matrices;
        const double *__restrict__ const charge_products;
        const double *__restrict__ const permutation_signs;
        const double *__restrict__ const permutation_matrices;

        double *__restrict__ const ax;
        double *__restrict__ const bx;
        double *__restrict__ const cx;
        double *__restrict__ const dx;

#ifdef ZSVM_SPHERICAL_ECG_CONTEXT_TIMING_ENABLED
        unsigned long long int matrix_element_calls;
        unsigned long long int matrix_element_time;
#endif // ZSVM_SPHERICAL_ECG_CONTEXT_TIMING_ENABLED

    private: // ============================================ FACTORY CONSTRUCTOR

        explicit SphericalECGContext(
                const std::vector<Particle> &particles,
                std::size_t num_permutations,
                int space_dimension,
                const Eigen::MatrixXd &inverse_mass_matrix,
                const Eigen::MatrixXd &pairwise_weight_vectors,
                const std::vector<double> &permutation_sign_vector,
                const std::vector<Eigen::MatrixXd> &permutation_matrix_vector);

    public: // ====================================================== DESTRUCTOR

        SphericalECGContext(const SphericalECGContext &) = delete;

        SphericalECGContext &operator=(const SphericalECGContext &) = delete;

        ~SphericalECGContext();

    public: // =========================================== STATIC FACTORY METHOD

        static SphericalECGContext create(
                const std::vector<Particle> &particles, int space_dimension);

    public: // ======================================================= ACCESSORS

#ifdef ZSVM_SPHERICAL_ECG_CONTEXT_TIMING_ENABLED

        unsigned long long int get_matrix_element_calls();

        unsigned long long int get_matrix_element_time();

#endif// ZSVM_SPHERICAL_ECG_CONTEXT_TIMING_ENABLED

    public: // =================================== MATRIX ELEMENT HELPER METHODS

        void gaussian_parameter_matrix(
                const double *__restrict__ correlation_coefficients,
                double *__restrict__ result) const;

        void matrix_element_kernel(
                double &__restrict__ overlap_kernel,
                double &__restrict__ hamiltonian_kernel);

    public: // ========================================== MATRIX ELEMENT METHODS

        void compute_matrix_elements(
                double &__restrict__ overlap_element,
                double &__restrict__ hamiltonian_element,
                const double *a, const double *b);

    public: // ========================================= RANDOM STATE GENERATION

        static std::mt19937_64 properly_seeded_random_engine();

        void random_correlation_coefficients(double *__restrict__ result) const;

    }; // class SphericalECGContext

} // namespace zsvm

#endif // ZSVM_SPHERICAL_ECG_CONTEXT_HPP_INCLUDED
