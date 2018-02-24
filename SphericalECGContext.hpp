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

        const Eigen::MatrixXd inverse_mass_matrix;
        const Eigen::MatrixXd pairwise_weight_vectors;
        const std::vector<Eigen::MatrixXd> pairwise_weight_matrices;
        const std::vector<double> pairwise_charge_products;
        const std::vector<double> permutation_signs;
        const std::vector<Eigen::MatrixXd> permutation_matrices;

    private: // ============================================ FACTORY CONSTRUCTOR

        explicit SphericalECGContext(
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

    public: // =================================== MATRIX ELEMENT HELPER METHODS

        Eigen::MatrixXd gaussian_parameter_matrix(
                const std::vector<double> &correlation_coefficients) const;

        void matrix_element_kernel(
                double &overlap_kernel, double &hamiltonian_kernel,
                const Eigen::MatrixXd &a, const Eigen::MatrixXd &b) const;

    public: // ========================================== MATRIX ELEMENT METHODS

        void matrix_elements(
                double &overlap_matrix_element,
                double &hamiltonian_matrix_element,
                const Eigen::MatrixXd &a, const Eigen::MatrixXd &b) const;

    public: // ========================================= RANDOM STATE GENERATION

        static std::mt19937_64 properly_seeded_random_engine();

        std::vector<double> random_correlation_coefficients();

    }; // class SphericalECGContext

} // namespace zsvm

#endif // ZSVM_SPHERICAL_ECG_CONTEXT_HPP_INCLUDED
