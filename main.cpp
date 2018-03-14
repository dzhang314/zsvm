// C++ standard library headers
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>

// Project-specific headers
#include "Particle.hpp"
#include "SphericalECGContext.hpp"
#include "RealVariationalSolver.hpp"
#include "AmoebaOptimizer.hpp"
#include "ScriptParser.hpp"

#define PRINT_EXECUTION_TIME(...) do { \
    ::std::chrono::high_resolution_clock::time_point __start_time = \
            ::std::chrono::high_resolution_clock::now(); \
    __VA_ARGS__; \
    ::std::chrono::high_resolution_clock::time_point __stop_time = \
            ::std::chrono::high_resolution_clock::now(); \
    ::std::cout << "Time elapsed (nanoseconds): " \
                << ::std::chrono::duration_cast<::std::chrono::nanoseconds>( \
                           __stop_time - __start_time).count() << std::endl; \
} while (0)


namespace zsvm {

    class SphericalECGVariationalOptimizer {

    public: // ================================================ MEMBER VARIABLES
        // TODO: public for debugging

        const std::size_t num_particles;
        const std::size_t num_pairs;
        SphericalECGContext context;
        RealVariationalSolver solver;
        std::vector<std::vector<double>> basis;
        std::vector<std::vector<double>> basis_matrices;

    public: // ===================================================== CONSTRUCTOR

        explicit SphericalECGVariationalOptimizer(
                const std::vector<Particle> &particles, int space_dimension)
                : num_particles(particles.size()),
                  num_pairs(num_particles * (num_particles - 1) / 2),
                  context(SphericalECGContext::create(
                          particles, space_dimension)) {}

    public: // =================================================================

        double get_ground_state_energy() {
            return solver.get_eigenvalue(0);
        }

    public: // =================================================================

        double augmented_ground_state_energy(
                const double *__restrict__ new_basis_matrix) {
            if (solver.empty()) {
                double overlap_element, hamiltonian_element;
                context.compute_matrix_elements(
                        overlap_element, hamiltonian_element,
                        new_basis_matrix, new_basis_matrix);
                return hamiltonian_element / overlap_element;
            } else {
                const std::size_t basis_size = basis_matrices.size();
                std::vector<double> new_overlap_column(basis_size + 1);
                std::vector<double> new_hamiltonian_column(basis_size + 1);
                for (std::size_t i = 0; i < basis_size; ++i) {
                    context.compute_matrix_elements(
                            new_overlap_column[i], new_hamiltonian_column[i],
                            basis_matrices[i].data(), new_basis_matrix);
                }
                context.compute_matrix_elements(
                        new_overlap_column[basis_size],
                        new_hamiltonian_column[basis_size],
                        new_basis_matrix, new_basis_matrix);
                return solver.minimum_augmented_eigenvalue(
                        new_overlap_column.data(),
                        new_hamiltonian_column.data());
            }
        }

        void expand_random(std::size_t num_trials) {
            std::vector<double> new_basis_element(num_pairs);
            std::vector<double> new_basis_matrix(num_pairs);
            std::vector<double> best_basis_element(num_pairs);
            std::vector<double> best_basis_matrix(num_pairs);
            double best_energy = solver.empty()
                                 ? std::numeric_limits<double>::max()
                                 : solver.get_eigenvalue(0);
            for (std::size_t trial = 0; trial < num_trials; ++trial) {
                context.random_correlation_coefficients(
                        new_basis_element.data());
                context.gaussian_parameter_matrix(
                        new_basis_element.data(), new_basis_matrix.data());
                const double new_energy =
                        augmented_ground_state_energy(new_basis_matrix.data());
                if (new_energy < best_energy) {
                    best_basis_element = new_basis_element;
                    best_basis_matrix = new_basis_matrix;
                    best_energy = new_energy;
                }
            }
            basis.push_back(best_basis_element);
            basis_matrices.push_back(best_basis_matrix);
            recompute_solver_matrices();
        }

        void expand_amoeba(std::size_t num_trials,
                           double initial_step_size,
                           std::size_t max_steps) {
            std::vector<double> new_basis_element(num_pairs);
            std::vector<double> best_basis_element(num_pairs);
            double best_energy = solver.empty()
                                 ? std::numeric_limits<double>::max()
                                 : solver.get_eigenvalue(0);
            for (std::size_t trial = 0; trial < num_trials; ++trial) {
                context.random_correlation_coefficients(
                        new_basis_element.data());
                const double new_energy = refine_amoeba(
                        new_basis_element.data(), initial_step_size, max_steps);
                if (new_energy < best_energy) {
                    std::cout << "Improved energy to "
                              << new_energy << std::endl;
                    best_basis_element = new_basis_element;
                    best_energy = new_energy;
                }
            }
            basis.push_back(best_basis_element);
            std::vector<double> best_basis_matrix(num_pairs);
            context.gaussian_parameter_matrix(
                    best_basis_element.data(), best_basis_matrix.data());
            basis_matrices.push_back(best_basis_matrix);
            recompute_solver_matrices();
        }

    public: // =================================================================

        double refine_amoeba(double *__restrict__ basis_element,
                             double initial_step_size,
                             std::size_t max_steps) {
            std::vector<double> basis_matrix(num_pairs);
            dznl::AmoebaOptimizer amoeba(
                    basis_element, num_pairs, initial_step_size,
                    [&](const double *b) {
                        context.gaussian_parameter_matrix(
                                b, basis_matrix.data());
                        return augmented_ground_state_energy(
                                basis_matrix.data());
                    });
            for (std::size_t step = 0; step < max_steps; ++step) {
                amoeba.step();
            }
            return amoeba.current_minimum(basis_element);
        }

    public: // =================================================================

        void recompute_basis_matrices() {
            basis_matrices.clear();
            for (const auto &basis_element : basis) {
                std::vector<double> basis_matrix(num_pairs);
                context.gaussian_parameter_matrix(
                        basis_element.data(), basis_matrix.data());
                basis_matrices.push_back(basis_matrix);
            }
        }

        void recompute_solver_matrices() {
            solver.set_basis_size_destructive(basis_matrices.size());
            for (std::size_t i = 0; i < basis_matrices.size(); ++i) {
                for (std::size_t j = 0; j < basis_matrices.size(); ++j) {
                    context.compute_matrix_elements(
                            solver.overlap_matrix_element(i, j),
                            solver.hamiltonian_matrix_element(i, j),
                            basis_matrices[i].data(),
                            basis_matrices[j].data());
                }
            }
        }

    }; // class SphericalECGVariationalOptimizer

} // namespace zsvm


int main() {
    std::cout << std::scientific;
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);

    const zsvm::Particle electron_up = {0, 1.0, -1.0, zsvm::Spin::UP};
    const zsvm::Particle electron_down = {0, 1.0, -1.0, zsvm::Spin::DOWN};
//    const zsvm::Particle positron_up = {1, 1.0, +1.0, zsvm::Spin::UP};
//    const zsvm::Particle positron_down = {1, 1.0, +1.0, zsvm::Spin::DOWN};
    const zsvm::Particle beryllium_nucleus =
            {2, 16538.028978017737, +4.0, zsvm::Spin::UP};
    const std::vector<zsvm::Particle> particles = {
            electron_up, electron_down, electron_up, electron_down,
            beryllium_nucleus};
    zsvm::SphericalECGVariationalOptimizer optimizer(particles, 3);

    for (std::size_t basis_size = 0; basis_size < 100; ++basis_size) {
        PRINT_EXECUTION_TIME(
                optimizer.expand_amoeba(10, 0.5, 200);
                optimizer.recompute_solver_matrices();
                std::cout << optimizer.get_ground_state_energy() << std::endl;
                std::ostringstream basis_output_file_name;
                basis_output_file_name << "basis-";
                basis_output_file_name << std::setw(8) << std::setfill('0');
                basis_output_file_name << basis_size + 1 << ".tsv";
                std::ofstream basis_output_file(basis_output_file_name.str());
                basis_output_file << std::scientific;
                basis_output_file << std::setprecision(
                        std::numeric_limits<double>::max_digits10);
                for (const auto &basis_element : optimizer.basis) {
                    for (std::size_t i = 0; i < basis_element.size(); ++i) {
                        if (i > 0) { basis_output_file << '\t'; }
                        basis_output_file << basis_element[i];
                    }
                    basis_output_file << std::endl;
                });
    }

#ifdef ZSVM_SPHERICAL_ECG_CONTEXT_TIMING_ENABLED
    std::cout << "Matrix element calls:       "
              << optimizer.context.get_matrix_element_calls() << std::endl;
    std::cout << "Matrix element time:        "
              << optimizer.context.get_matrix_element_time() << std::endl;
#endif // ZSVM_SPHERICAL_ECG_CONTEXT_TIMING_ENABLED

    return 0;
}
