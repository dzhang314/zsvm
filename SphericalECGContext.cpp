#include "SphericalECGContext.hpp"

// Project-specific headers
#include "PackedLinearAlgebra.hpp" // for packed_determinant_inverse et al.


const zsvm::packed_determinant_inverse_function_t
        zsvm::SphericalECGContext::PACKED_DETERMINANT_INVERSE[] = {
        packed_determinant_inverse_1,
        packed_determinant_inverse_2,
        packed_determinant_inverse_3,
        packed_determinant_inverse_4,
        packed_determinant_inverse_5,
};


const zsvm::packed_kinetic_trace_function_t
        zsvm::SphericalECGContext::PACKED_KINETIC_TRACE[] = {
        packed_kinetic_trace_1,
        packed_kinetic_trace_2,
        packed_kinetic_trace_3,
        packed_kinetic_trace_4,
        packed_kinetic_trace_5,
};


const zsvm::packed_quadratic_form_function_t
        zsvm::SphericalECGContext::PACKED_QUADRATIC_FORM[] = {
        packed_quadratic_form_1,
        packed_quadratic_form_2,
        packed_quadratic_form_3,
        packed_quadratic_form_4,
        packed_quadratic_form_5,
};


const zsvm::packed_permutation_conjugate_function_t
        zsvm::SphericalECGContext::PACKED_PERMUTATION_CONJUGATE[] = {
        packed_permutation_conjugate_1,
        packed_permutation_conjugate_2,
        packed_permutation_conjugate_3,
        packed_permutation_conjugate_4,
        packed_permutation_conjugate_5,
};


//// C++ standard library headers
//#include <chrono> // for std::chrono::system_clock
//#include <cmath> // for std::exp, std::tgamma
//#include <stdexcept> // for std::invalid_argument
//
//// Eigen linear algebra library headers
//#include <Eigen/LU> // for inverse(), determinant()
//
//#include "Permutation.hpp" // for dznl::invariant_permutations et al.
//#include "JacobiCoordinates.hpp" // for jaco::pairwise_weights et al.
//
//zsvm::SphericalECGContext::SphericalECGContext(
//        const std::vector<Particle> &particles,
//        std::size_t num_permutations,
//        int space_dimension,
//        const Eigen::MatrixXd &inverse_mass_matrix,
//        const Eigen::MatrixXd &pairwise_weight_vectors,
//        const std::vector<double> &permutation_sign_vector,
//        const std::vector<Eigen::MatrixXd> &permutation_matrix_vector)
//        : num_particles(particles.size()),
//          num_pairs(num_particles * (num_particles - 1) / 2),
//          num_permutations(num_permutations),
//          matrix_size((num_particles - 1) * (num_particles - 1)),
//          space_dimension(space_dimension),
//          dimension_factor(std::tgamma(0.5 * (space_dimension - 1.0)) /
//                           std::tgamma(0.5 * space_dimension)),
//          kinetic_factor((0.5 * space_dimension) / dimension_factor),
//          inverse_masses(new double[num_particles - 1]),
//          weight_vectors(new double[num_pairs * (num_particles - 1)]),
//          weight_matrices(new double[num_pairs * num_pairs]),
//          charge_products(new double[num_pairs]),
//          permutation_signs(new double[num_permutations]),
//          permutation_matrices(new double[num_permutations * matrix_size]),
//          ax(new double[num_pairs]),
//          bx(new double[num_pairs]),
//          cx(new double[num_pairs]),
//          dx(new double[num_pairs]),
//          packed_determinant_inverse(
//                  PACKED_DETERMINANT_INVERSE[num_particles - 2]),
//          packed_kinetic_trace(
//                  PACKED_KINETIC_TRACE[num_particles - 2]),
//          packed_quadratic_form(
//                  PACKED_QUADRATIC_FORM[num_particles - 2]),
//          packed_permutation_conjugate(
//                  PACKED_PERMUTATION_CONJUGATE[num_particles - 2])
//#ifdef ZSVM_SPHERICAL_ECG_CONTEXT_TIMING_ENABLED
//        , matrix_element_calls(0), matrix_element_time(0)
//#endif // ZSVM_SPHERICAL_ECG_CONTEXT_TIMING_ENABLED
//{
//    auto inverse_mass_pointer = const_cast<double *>(inverse_masses);
//    for (std::size_t i = 0; i < num_particles - 1; ++i) {
//        inverse_mass_pointer[i] = inverse_mass_matrix(i, i);
//    }
//    auto weight_vector_pointer = const_cast<double *>(weight_vectors);
//    for (std::size_t p = 0, k = 0; p < num_pairs; ++p) {
//        for (std::size_t i = 0; i < num_particles - 1; ++i, ++k) {
//            weight_vector_pointer[k] = pairwise_weight_vectors(i, p);
//        }
//    }
//    auto weight_matrix_pointer = const_cast<double *>(weight_matrices);
//    for (std::size_t p = 0, k = 0; p < num_pairs; ++p) {
//        for (std::size_t i = 0; i < num_particles - 1; ++i) {
//            for (std::size_t j = 0; j <= i; ++j, ++k) {
//                weight_matrix_pointer[k] = pairwise_weight_vectors(i, p) *
//                                           pairwise_weight_vectors(j, p);
//            }
//        }
//    }
//    auto charge_product_pointer = const_cast<double *>(charge_products);
//    for (std::size_t i = 0, k = 0; i < num_particles - 1; ++i) {
//        for (std::size_t j = i + 1; j < num_particles; ++j, ++k) {
//            charge_product_pointer[k] = particles[i].charge *
//                                        particles[j].charge;
//        }
//    }
//    auto permutation_sign_pointer = const_cast<double *>(permutation_signs);
//    for (std::size_t p = 0; p < num_permutations; ++p) {
//        permutation_sign_pointer[p] = permutation_sign_vector[p];
//    }
//    auto permutation_matrix_pointer =
//            const_cast<double *>(permutation_matrices);
//    for (std::size_t p = 0, k = 0; p < num_permutations; ++p) {
//        for (std::size_t i = 0; i < num_particles - 1; ++i) {
//            for (std::size_t j = 0; j < num_particles - 1; ++j, ++k) {
//                permutation_matrix_pointer[k] =
//                        permutation_matrix_vector[p](j, i);
//            }
//        }
//    }
//}
//
//
//zsvm::SphericalECGContext::~SphericalECGContext() {
//    delete[] inverse_masses;
//    delete[] weight_vectors;
//    delete[] weight_matrices;
//    delete[] charge_products;
//    delete[] permutation_signs;
//    delete[] permutation_matrices;
//    delete[] ax;
//    delete[] bx;
//    delete[] cx;
//    delete[] dx;
//}
//
//
//zsvm::SphericalECGContext zsvm::SphericalECGContext::create(
//        const std::vector<Particle> &particles, int space_dimension) {
//    const std::size_t num_particles = particles.size();
//    if (num_particles < 2) {
//        throw std::invalid_argument(
//                "Attempted to construct SVMSolver "
//                        "with fewer than 2 particles");
//    }
//    std::vector<int> particle_types(num_particles);
//    std::vector<double> masses(num_particles);
//    std::vector<Spin> spins(num_particles);
//    for (std::size_t i = 0; i < num_particles; ++i) {
//        particle_types[i] = particles[i].type;
//        masses[i] = particles[i].mass;
//        spins[i] = particles[i].spin;
//    }
//    const std::vector<std::vector<std::size_t>> allowed_permutations =
//            dznl::invariant_permutations(particle_types);
//    std::vector<double> permutation_signs;
//    for (const auto &permutation : allowed_permutations) {
//        const std::size_t signature =
//                dznl::count_changes(spins, permutation) / 2 +
//                dznl::count_inversions(permutation);
//        permutation_signs.push_back((signature % 2 == 0) ? +1.0 : -1.0);
//    }
//    return SphericalECGContext(
//            particles,
//            allowed_permutations.size(),
//            space_dimension,
//            jaco::reduced_inverse_mass_matrix(masses),
//            jaco::pairwise_weights(masses),
//            permutation_signs,
//            jaco::permutation_matrices(masses, allowed_permutations));
//}
//
//
//#ifdef ZSVM_SPHERICAL_ECG_CONTEXT_TIMING_ENABLED
//
//unsigned long long int zsvm::SphericalECGContext::get_matrix_element_calls() {
//    return matrix_element_calls;
//}
//
//#endif // ZSVM_SPHERICAL_ECG_CONTEXT_TIMING_ENABLED
//
//
//#ifdef ZSVM_SPHERICAL_ECG_CONTEXT_TIMING_ENABLED
//
//unsigned long long int zsvm::SphericalECGContext::get_matrix_element_time() {
//    return matrix_element_time;
//}
//
//#endif // ZSVM_SPHERICAL_ECG_CONTEXT_TIMING_ENABLED
//
//
//void zsvm::SphericalECGContext::gaussian_parameter_matrix(
//        const double *__restrict__ correlation_coefficients,
//        double *__restrict__ result) const {
//    for (std::size_t p = 0; p < num_pairs; ++p) { result[p] = 0.0; }
//    for (std::size_t p = 0, k = 0; p < num_pairs; ++p) {
//        const double c = std::exp(correlation_coefficients[p]);
//        for (std::size_t q = 0; q < num_pairs; ++q, ++k) {
//            result[q] += c * weight_matrices[k];
//        }
//    }
//}
//
//
//static inline double half_inverse_pow(double x, int n) {
//    switch (n) {
//        case 0:
//            return 1.0;
//        case 1:
//            return 1.0 / std::sqrt(x);
//        case 2:
//            return 1.0 / x;
//        case 3:
//            return 1.0 / (x * std::sqrt(x));
//        case 4:
//            return 1.0 / (x * x);
//        case 5:
//            return 1.0 / (x * x * std::sqrt(x));
//        case 6:
//            return 1.0 / (x * x * x);
//        case 7:
//            return 1.0 / (x * x * x * std::sqrt(x));
//        default:
//            return 1.0 / std::sqrt(std::pow(x, n));
//    }
//}
//
//
//void zsvm::SphericalECGContext::matrix_element_kernel(
//        double &__restrict__ overlap_kernel,
//        double &__restrict__ hamiltonian_kernel) {
//    for (std::size_t i = 0; i < num_pairs; ++i) { cx[i] = ax[i] + bx[i]; }
//    overlap_kernel = half_inverse_pow(
//            packed_determinant_inverse(cx, dx), space_dimension);
//    hamiltonian_kernel = kinetic_factor * packed_kinetic_trace(
//            ax, bx, dx, inverse_masses);
//    for (std::size_t k = 0; k < num_pairs; ++k) {
//        const double alpha = packed_quadratic_form(
//                dx, weight_vectors + (num_particles - 1) * k);
//        hamiltonian_kernel += charge_products[k] / std::sqrt(2.0 * alpha);
//    }
//    hamiltonian_kernel *= dimension_factor * overlap_kernel;
//}
//
//
//void zsvm::SphericalECGContext::compute_matrix_elements(
//        double &__restrict__ overlap_element,
//        double &__restrict__ hamiltonian_element,
//        const double *a, const double *b) {
//#ifdef ZSVM_SPHERICAL_ECG_CONTEXT_TIMING_ENABLED
//    auto start = std::chrono::high_resolution_clock::now();
//#endif // ZSVM_SPHERICAL_ECG_CONTEXT_TIMING_ENABLED
//    overlap_element = hamiltonian_element = 0.0;
//    double overlap_kernel, hamiltonian_kernel;
//    for (std::size_t i = 0; i < num_permutations; ++i) {
//        for (std::size_t j = 0; j < num_permutations; ++j) {
//            const double sign = permutation_signs[i] * permutation_signs[j];
//            packed_permutation_conjugate(
//                    a, permutation_matrices + i * matrix_size, ax);
//            packed_permutation_conjugate(
//                    b, permutation_matrices + j * matrix_size, bx);
//            matrix_element_kernel(overlap_kernel, hamiltonian_kernel);
//            overlap_element += sign * overlap_kernel;
//            hamiltonian_element += sign * hamiltonian_kernel;
//        }
//    }
//#ifdef ZSVM_SPHERICAL_ECG_CONTEXT_TIMING_ENABLED
//    auto stop = std::chrono::high_resolution_clock::now();
//    matrix_element_time += std::chrono::duration_cast<
//            std::chrono::nanoseconds>(stop - start).count();
//    ++matrix_element_calls;
//#endif // ZSVM_SPHERICAL_ECG_CONTEXT_TIMING_ENABLED
//}
//
//
//std::mt19937_64 zsvm::SphericalECGContext::properly_seeded_random_engine() {
//    std::array<std::mt19937_64::result_type,
//            std::mt19937_64::state_size> random_data;
//    std::random_device source;
//    std::generate(random_data.begin(), random_data.end(), std::ref(source));
//    random_data[0] = static_cast<std::mt19937_64::result_type>(
//            std::chrono::duration_cast<std::chrono::milliseconds>(
//                    std::chrono::system_clock::now().time_since_epoch()
//            ).count());
//    std::seed_seq seeds(random_data.begin(), random_data.end());
//    return std::mt19937_64(seeds);
//}
//
//
//void zsvm::SphericalECGContext::random_correlation_coefficients(
//        double *__restrict__ result) const {
//    static std::mt19937_64 random_engine = properly_seeded_random_engine();
//    static std::normal_distribution<double> correlation_distribution(0.0, 3.0);
//    for (std::size_t i = 0; i < num_pairs; ++i) {
//        result[i] = correlation_distribution(random_engine);
//    }
//}
