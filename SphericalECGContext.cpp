#include "SphericalECGContext.hpp"

// C++ standard library headers
#include <chrono> // for std::chrono::system_clock
#include <cmath> // for std::exp, std::tgamma
#include <stdexcept> // for std::invalid_argument

// Eigen linear algebra library headers
#include <Eigen/LU>

// Project-specific headers
#include "Permutation.hpp"
#include "JacobiCoordinates.hpp"

zsvm::SphericalECGContext::SphericalECGContext(
        std::size_t num_pairs,
        std::size_t num_permutations,
        int space_dimension,
        double dimension_factor,
        const Eigen::MatrixXd &inverse_mass_matrix,
        const Eigen::MatrixXd &pairwise_weight_vectors,
        const std::vector<Eigen::MatrixXd> &pairwise_weight_matrices,
        const std::vector<double> &pairwise_charge_products,
        const std::vector<double> &permutation_signs,
        const std::vector<Eigen::MatrixXd> &permutation_matrices)
        : num_pairs(num_pairs),
          num_permutations(num_permutations),
          space_dimension(space_dimension),
          dimension_factor(dimension_factor),
          inverse_mass_matrix(inverse_mass_matrix),
          pairwise_weight_vectors(pairwise_weight_vectors),
          pairwise_weight_matrices(pairwise_weight_matrices),
          pairwise_charge_products(pairwise_charge_products),
          permutation_signs(permutation_signs),
          permutation_matrices(permutation_matrices) {}

zsvm::SphericalECGContext zsvm::SphericalECGContext::create(
        const std::vector<Particle> &particles, int space_dimension) {

    const std::size_t num_particles = particles.size();
    if (num_particles < 2) {
        throw std::invalid_argument(
                "Attempted to construct SVMSolver "
                        "with fewer than 2 particles");
    }

    std::vector<int> particle_types(num_particles);
    std::vector<double> masses(num_particles);
    std::vector<double> charges(num_particles);
    std::vector<Spin> spins(num_particles);
    for (std::size_t i = 0; i < num_particles; ++i) {
        particle_types[i] = particles[i].type;
        masses[i] = particles[i].mass;
        charges[i] = particles[i].charge;
        spins[i] = particles[i].spin;
    }

    const std::size_t num_pairs =
            num_particles * (num_particles - 1) / 2;
    const Eigen::MatrixXd pairwise_weight_vectors =
            jaco::pairwise_weights(masses);
    std::vector<Eigen::MatrixXd> pairwise_weight_matrices;
    for (std::size_t k = 0; k < num_pairs; ++k) {
        pairwise_weight_matrices.push_back(
                pairwise_weight_vectors.col(k) *
                pairwise_weight_vectors.col(k).transpose());
    }

    std::vector<double> pairwise_charge_products(num_pairs);
    for (std::size_t i = 0, k = 0; i < num_particles - 1; ++i) {
        for (std::size_t j = i + 1; j < num_particles; ++j, ++k) {
            pairwise_charge_products[k] = charges[i] * charges[j];
        }
    }

    const std::vector<std::vector<std::size_t>> allowed_permutations =
            dznl::invariant_permutations(particle_types);
    std::vector<double> permutation_signs;
    for (const auto &permutation : allowed_permutations) {
        const std::size_t signature =
                dznl::count_changes(spins, permutation) / 2 +
                dznl::count_inversions(permutation);
        permutation_signs.push_back((signature % 2 == 0)
                                    ? +1.0 : -1.0);
    }

    return SphericalECGContext(
            num_pairs,
            allowed_permutations.size(),
            space_dimension,
            std::tgamma(0.5 * (space_dimension - 1.0)) /
            std::tgamma(0.5 * space_dimension),
            jaco::reduced_inverse_mass_matrix(masses),
            pairwise_weight_vectors,
            pairwise_weight_matrices,
            pairwise_charge_products,
            permutation_signs,
            jaco::permutation_matrices(masses, allowed_permutations));
}

Eigen::MatrixXd zsvm::SphericalECGContext::gaussian_parameter_matrix(
        const std::vector<double> &correlation_coefficients) const {
    if (correlation_coefficients.size() != num_pairs) {
        throw std::invalid_argument(
                "SphericalECGContext::gaussian_parameter_matrix "
                        "received vector with incorrect number of "
                        "correlation coefficients");
    }
    Eigen::MatrixXd a =
            correlation_coefficients[0] * pairwise_weight_matrices[0];
    for (std::size_t k = 1; k < num_pairs; ++k) {
        a += correlation_coefficients[k] * pairwise_weight_matrices[k];
    }
    return a;
}

void zsvm::SphericalECGContext::matrix_element_kernel(
        double &overlap_kernel, double &hamiltonian_kernel,
        const Eigen::MatrixXd &a, const Eigen::MatrixXd &b) const {
    const Eigen::MatrixXd c = a + b;
    const Eigen::MatrixXd c_inv = c.inverse();
    overlap_kernel = std::sqrt(
            std::pow(c.determinant(), -space_dimension));
    hamiltonian_kernel = (0.5 * space_dimension) / dimension_factor *
                         (inverse_mass_matrix * a * c_inv * b).trace();
    for (std::size_t k = 0; k < num_pairs; ++k) {
        hamiltonian_kernel += pairwise_charge_products[k] / std::sqrt(
                2.0 * pairwise_weight_vectors.col(k).dot(
                        c_inv * pairwise_weight_vectors.col(k)));
    }
    hamiltonian_kernel *= dimension_factor * overlap_kernel;
}

void zsvm::SphericalECGContext::matrix_elements(
        double &overlap_matrix_element,
        double &hamiltonian_matrix_element,
        const Eigen::MatrixXd &a, const Eigen::MatrixXd &b) const {
    overlap_matrix_element = hamiltonian_matrix_element = 0.0;
    double overlap_kernel, hamiltonian_kernel;
    for (std::size_t i = 0; i < num_permutations; ++i) {
        for (std::size_t j = 0; j < num_permutations; ++j) {
            const double sign = permutation_signs[i] *
                                permutation_signs[j];
            matrix_element_kernel(overlap_kernel, hamiltonian_kernel,
                                  permutation_matrices[i].transpose() *
                                  a * permutation_matrices[i],
                                  permutation_matrices[j].transpose() *
                                  b * permutation_matrices[j]);
            overlap_matrix_element += sign * overlap_kernel;
            hamiltonian_matrix_element += sign * hamiltonian_kernel;
        }
    }
}

std::mt19937_64 zsvm::SphericalECGContext::properly_seeded_random_engine() {
    std::array<std::mt19937_64::result_type,
            std::mt19937_64::state_size> random_data;
    std::random_device source;
    std::generate(random_data.begin(), random_data.end(),
                  std::ref(source));
    random_data[0] = static_cast<std::mt19937_64::result_type>(
            std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::system_clock::now().time_since_epoch()
            ).count());
    std::seed_seq seeds(random_data.begin(), random_data.end());
    return std::mt19937_64(seeds);
}

std::vector<double>
zsvm::SphericalECGContext::random_correlation_coefficients() {
    static std::mt19937_64 random_engine =
            properly_seeded_random_engine();
    static std::normal_distribution<double>
            correlation_distribution(0.0, 8.0);
    std::vector<double> result;
    for (std::size_t i = 0; i < num_pairs; ++i) {
        result.push_back(std::exp(
                correlation_distribution(random_engine)));
    }
    return result;
}
