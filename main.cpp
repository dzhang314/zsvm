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
    ::std::cout << "Time elapsed (microseconds): " \
                << ::std::chrono::duration_cast<::std::chrono::microseconds>( \
                           __stop_time - __start_time).count() << std::endl; \
} while (0)


int main() {
    const zsvm::Particle electron_up = {0, 1.0, -1.0, zsvm::Spin::UP};
    const zsvm::Particle electron_down = {0, 1.0, -1.0, zsvm::Spin::DOWN};
//    const zsvm::Particle positron_up = {1, 1.0, +1.0, zsvm::Spin::UP};
//    const zsvm::Particle positron_down = {1, 1.0, +1.0, zsvm::Spin::DOWN};
    const zsvm::Particle beryllium_nucleus = {2, 16538.028978017737,
                                              +4.0, zsvm::Spin::UP};
    const std::vector<zsvm::Particle> particles = {
            electron_up, electron_down, electron_up, electron_down,
            beryllium_nucleus};
    const zsvm::SphericalECGContext context =
            zsvm::SphericalECGContext::create(particles, 3);

    std::size_t basis_size = 2;
    std::vector<std::vector<double>> basis;
    for (std::size_t i = 0; i < basis_size; ++i) {
        basis.push_back(context.random_correlation_coefficients());
    }
    std::vector<Eigen::MatrixXd> basis_matrices;
    for (const auto &correlation_coefficients : basis) {
        basis_matrices.push_back(
                context.gaussian_parameter_matrix(correlation_coefficients));
    }
    zsvm::RealVariationalSolver solver;
    solver.set_basis_size_destructive(basis_size);
    for (std::size_t i = 0; i < basis_size; ++i) {
        for (std::size_t j = 0; j < basis_size; ++j) {
            double overlap_matrix_element, hamiltonian_matrix_element;
            context.compute_matrix_elements(
                    overlap_matrix_element, hamiltonian_matrix_element,
                    basis_matrices[i], basis_matrices[j]);
            solver.set_overlap_matrix_element(
                    i, j, overlap_matrix_element);
            solver.set_hamiltonian_matrix_element(
                    i, j, hamiltonian_matrix_element);
        }
    }

    for (std::size_t basis_index = 0; basis_index < 100; ++basis_index) {
        std::cout << "Initial energy: " << solver.get_eigenvalue(0) << std::endl;
        std::vector<double> best_basis_element;
        Eigen::MatrixXd best_basis_matrix;
        double best_energy = solver.get_eigenvalue(0);
        for (std::size_t trial_index = 0; trial_index < 100; ++trial_index) {
            const std::vector<double> new_basis_element =
                    context.random_correlation_coefficients();
            const Eigen::MatrixXd new_basis_matrix =
                    context.gaussian_parameter_matrix(new_basis_element);
            std::vector<double> new_overlap_column, new_hamiltonian_column;
            for (std::size_t i = 0; i < basis_size; ++i) {
                double overlap_matrix_element, hamiltonian_matrix_element;
                context.compute_matrix_elements(
                        overlap_matrix_element, hamiltonian_matrix_element,
                        basis_matrices[i], new_basis_matrix);
                new_overlap_column.push_back(overlap_matrix_element);
                new_hamiltonian_column.push_back(hamiltonian_matrix_element);
            }
            {
                double overlap_matrix_element, hamiltonian_matrix_element;
                context.compute_matrix_elements(
                        overlap_matrix_element, hamiltonian_matrix_element,
                        new_basis_matrix, new_basis_matrix);
                new_overlap_column.push_back(overlap_matrix_element);
                new_hamiltonian_column.push_back(hamiltonian_matrix_element);
            }
            const double new_energy = solver.minimum_augmented_eigenvalue(
                    &new_overlap_column[0], &new_hamiltonian_column[0]);
            if (new_energy < best_energy) {
                std::cout << "Improved energy to: " << new_energy << std::endl;
                best_basis_element = new_basis_element;
                best_basis_matrix = new_basis_matrix;
                best_energy = new_energy;
            }
        }
        basis.push_back(best_basis_element);
        basis_matrices.push_back(best_basis_matrix);
        solver.set_basis_size_destructive(++basis_size);
        for (std::size_t i = 0; i < basis_size; ++i) {
            for (std::size_t j = 0; j < basis_size; ++j) {
                double overlap_matrix_element, hamiltonian_matrix_element;
                context.compute_matrix_elements(
                        overlap_matrix_element, hamiltonian_matrix_element,
                        basis_matrices[i], basis_matrices[j]);
                solver.set_overlap_matrix_element(
                        i, j, overlap_matrix_element);
                solver.set_hamiltonian_matrix_element(
                        i, j, hamiltonian_matrix_element);
            }
        }
    }
    return 0;
}
