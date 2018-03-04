// C++ standard library headers
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

// Project-specific headers
#include "Particle.hpp"
#include "SphericalECGContext.hpp"
#include "RealVariationalSolver.hpp"


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

        SphericalECGContext context;
        RealVariationalSolver solver;
        std::vector<std::vector<double>> basis;
        std::vector<std::vector<double>> basis_matrices;

    public: // ===================================================== CONSTRUCTOR

        explicit SphericalECGVariationalOptimizer(
                const std::vector<Particle> &particles, int space_dimension)
                : context(SphericalECGContext::create(particles,
                                                      space_dimension)) {}

    public: // =================================================================

        double get_ground_state_energy() {
            return solver.get_eigenvalue(0);
        }

    public: // =================================================================

        double augmented_ground_state_energy(
                const std::vector<double> &basis_matrix) {
            if (solver.empty()) {
                double overlap_element, hamiltonian_element;
                context.compute_matrix_elements(
                        overlap_element, hamiltonian_element,
                        basis_matrix.data(), basis_matrix.data());
                return hamiltonian_element / overlap_element;
            } else {
                std::vector<double> new_overlap_column, new_hamiltonian_column;
                for (std::size_t i = 0; i < basis_matrices.size(); ++i) {
                    double overlap_element, hamiltonian_element;
                    context.compute_matrix_elements(
                            overlap_element, hamiltonian_element,
                            basis_matrices[i].data(), basis_matrix.data());
                    new_overlap_column.push_back(overlap_element);
                    new_hamiltonian_column.push_back(hamiltonian_element);
                }
                double overlap_matrix_element, hamiltonian_matrix_element;
                context.compute_matrix_elements(
                        overlap_matrix_element, hamiltonian_matrix_element,
                        basis_matrix.data(), basis_matrix.data());
                new_overlap_column.push_back(overlap_matrix_element);
                new_hamiltonian_column.push_back(hamiltonian_matrix_element);
                return solver.minimum_augmented_eigenvalue(
                        &new_overlap_column[0], &new_hamiltonian_column[0]);
            }
        }

        void expand_stochastic(std::size_t num_trials) {
            std::vector<double> best_basis_element;
            std::vector<double> best_basis_matrix;
            double best_energy = solver.empty()
                                 ? std::numeric_limits<double>::max()
                                 : solver.get_eigenvalue(0);
            for (std::size_t trial = 0; trial < num_trials; ++trial) {
                const std::vector<double> new_basis_element =
                        context.random_correlation_coefficients();
                const std::vector<double> new_basis_matrix =
                        context.gaussian_parameter_matrix(new_basis_element);
                const double new_energy =
                        augmented_ground_state_energy(new_basis_matrix);
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

    public: // =================================================================

        void recompute_basis_matrices() {
            basis_matrices.clear();
            for (const auto &basis_element : basis) {
                basis_matrices.push_back(
                        context.gaussian_parameter_matrix(basis_element));
            }
        }

        void recompute_solver_matrix_elements(std::size_t i, std::size_t j) {
            double overlap_element, hamiltonian_element;
            context.compute_matrix_elements(
                    overlap_element, hamiltonian_element,
                    basis_matrices[i].data(), basis_matrices[j].data());
            solver.set_overlap_matrix_element(i, j, overlap_element);
            solver.set_hamiltonian_matrix_element(i, j, hamiltonian_element);
        }

        void recompute_solver_matrices() {
            solver.set_basis_size_destructive(basis_matrices.size());
            for (std::size_t i = 0; i < basis_matrices.size(); ++i) {
                for (std::size_t j = 0; j < basis_matrices.size(); ++j) {
                    recompute_solver_matrix_elements(i, j);
                }
            }
        }

    }; // class SphericalECGVariationalOptimizer

} // namespace zsvm


int main() {
    const zsvm::Particle electron_up =
            {0, 1.0, -1.0, zsvm::Spin::UP};
    const zsvm::Particle electron_down =
            {0, 1.0, -1.0, zsvm::Spin::DOWN};
//    const zsvm::Particle positron_up = {1, 1.0, +1.0, zsvm::Spin::UP};
//    const zsvm::Particle positron_down = {1, 1.0, +1.0, zsvm::Spin::DOWN};
    const zsvm::Particle beryllium_nucleus =
            {2, 16538.028978017737, +4.0, zsvm::Spin::UP};
    const std::vector<zsvm::Particle> particles = {
            electron_up, electron_down, electron_up, electron_down,
            beryllium_nucleus};
    zsvm::SphericalECGVariationalOptimizer optimizer(particles, 3);

    PRINT_EXECUTION_TIME(
            for (std::size_t basis_size = 0; basis_size < 10; ++basis_size) {
                optimizer.expand_stochastic(200);
                optimizer.recompute_solver_matrices();
                std::cout << optimizer.get_ground_state_energy() << std::endl;
            });

#ifdef ZSVM_SPHERICAL_ECG_CONTEXT_TIMING_ENABLED
    std::cout << "Matrix element calls:       "
              << optimizer.context.get_matrix_element_calls() << std::endl;
    std::cout << "Matrix element time:        "
              << optimizer.context.get_matrix_element_time() << std::endl;
#endif // ZSVM_SPHERICAL_ECG_CONTEXT_TIMING_ENABLED

    return 0;
}
