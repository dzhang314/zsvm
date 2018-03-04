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

namespace zsvm {

    class SphericalECGContext {

    private: // =============================================== MEMBER VARIABLES

        const std::size_t num_pairs;
        const std::size_t num_permutations;
        const int space_dimension;
        const double dimension_factor;
        const double kinetic_factor;

        const Eigen::MatrixXd inverse_mass_matrix;
        const Eigen::MatrixXd pairwise_weight_vectors;
        const std::vector<Eigen::MatrixXd> pairwise_weight_matrices;
        const std::vector<double> pairwise_charge_products;
        const std::vector<double> permutation_signs;
        const std::vector<Eigen::MatrixXd> permutation_matrices;

        Eigen::MatrixXd workspace;
        Eigen::MatrixXd workspace_inverse;

        unsigned long long int matrix_element_calls;
        unsigned long long int matrix_element_time;

    private: // ============================================ FACTORY CONSTRUCTOR

        explicit SphericalECGContext(
                std::size_t num_particles,
                std::size_t num_pairs,
                std::size_t num_permutations,
                int space_dimension,
                double dimension_factor,
                const Eigen::MatrixXd &inverse_mass_matrix,
                const Eigen::MatrixXd &pairwise_weight_vectors,
                const std::vector<Eigen::MatrixXd> &pairwise_weight_matrices,
                const std::vector<double> &pairwise_charge_products,
                const std::vector<double> &permutation_signs,
                const std::vector<Eigen::MatrixXd> &permutation_matrices);

    public: // =========================================== STATIC FACTORY METHOD

        static SphericalECGContext create(
                const std::vector<Particle> &particles, int space_dimension);

    public: // ======================================================= ACCESSORS

        unsigned long long int get_matrix_element_calls();

        unsigned long long int get_matrix_element_time();

    public: // =================================== MATRIX ELEMENT HELPER METHODS

        Eigen::MatrixXd gaussian_parameter_matrix(
                const std::vector<double> &correlation_coefficients) const;

        void matrix_element_kernel(
                double &__restrict__ overlap_kernel,
                double &__restrict__ hamiltonian_kernel,
                const Eigen::MatrixXd &__restrict__ a,
                const Eigen::MatrixXd &__restrict__ b);

    public: // ========================================== MATRIX ELEMENT METHODS

        void compute_matrix_elements(
                double &overlap_matrix_element,
                double &hamiltonian_matrix_element,
                const Eigen::MatrixXd &a, const Eigen::MatrixXd &b);

    public: // ========================================= RANDOM STATE GENERATION

        static std::mt19937_64 properly_seeded_random_engine();

        std::vector<double> random_correlation_coefficients() const;

    }; // class SphericalECGContext

} // namespace zsvm

#endif // ZSVM_SPHERICAL_ECG_CONTEXT_HPP_INCLUDED
