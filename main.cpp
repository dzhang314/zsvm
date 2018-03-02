// C++ standard library headers
#include <chrono>
#include <iomanip>
#include <iostream>
#include <thread>
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
    std::cout << std::setprecision(16) << std::endl;

    zsvm::RealVariationalSolver test_solver;
    test_solver.set_basis_size_destructive(2);
    test_solver.set_overlap_matrix_element(0, 0, +2.0);
    test_solver.set_overlap_matrix_element(0, 1, +1.0);
    test_solver.set_overlap_matrix_element(1, 0, +1.0);
    test_solver.set_overlap_matrix_element(1, 1, +2.0);
    test_solver.set_hamiltonian_matrix_element(0, 0, -3.0);
    test_solver.set_hamiltonian_matrix_element(0, 1, +4.0);
    test_solver.set_hamiltonian_matrix_element(1, 0, +4.0);
    test_solver.set_hamiltonian_matrix_element(1, 1, -5.0);

    std::vector<double> new_overlap_column = {0.5, 0.25, 3.0};
    std::vector<double> new_hamiltonian_column = {0.0, 2.0, 3.0};

    std::cout << test_solver.meuds(&new_overlap_column[0],
                                   &new_hamiltonian_column[0]);


    zsvm::RealVariationalSolver test_solver_2;
    test_solver_2.set_basis_size_destructive(3);
    test_solver_2.set_overlap_matrix_element(0, 0, 2.0);
    test_solver_2.set_overlap_matrix_element(0, 1, 1.0);
    test_solver_2.set_overlap_matrix_element(0, 2, 0.5);
    test_solver_2.set_overlap_matrix_element(1, 0, 1.0);
    test_solver_2.set_overlap_matrix_element(1, 1, 2.0);
    test_solver_2.set_overlap_matrix_element(1, 2, 0.25);
    test_solver_2.set_overlap_matrix_element(2, 0, 0.5);
    test_solver_2.set_overlap_matrix_element(2, 1, 0.25);
    test_solver_2.set_overlap_matrix_element(2, 2, 3.0);
    test_solver_2.set_hamiltonian_matrix_element(0, 0, -3.0);
    test_solver_2.set_hamiltonian_matrix_element(0, 1, 4.0);
    test_solver_2.set_hamiltonian_matrix_element(0, 2, 0.0);
    test_solver_2.set_hamiltonian_matrix_element(1, 0, 4.0);
    test_solver_2.set_hamiltonian_matrix_element(1, 1, -5.0);
    test_solver_2.set_hamiltonian_matrix_element(1, 2, 2.0);
    test_solver_2.set_hamiltonian_matrix_element(2, 0, 0.0);
    test_solver_2.set_hamiltonian_matrix_element(2, 1, 2.0);
    test_solver_2.set_hamiltonian_matrix_element(2, 2, 3.0);
    std::cout << test_solver_2.get_eigenvalue(0) << std::endl;

//    const zsvm::Particle electron_up = {0, 1.0, -1.0, zsvm::Spin::UP};
//    const zsvm::Particle electron_down = {0, 1.0, -1.0, zsvm::Spin::DOWN};
////    const zsvm::Particle positron_up = {1, 1.0, +1.0, zsvm::Spin::UP};
////    const zsvm::Particle positron_down = {1, 1.0, +1.0, zsvm::Spin::DOWN};
//    const zsvm::Particle beryllium_nucleus = {2, 16538.028978017737,
//                                              +4.0, zsvm::Spin::UP};
//    const std::vector<zsvm::Particle> particles = {
//            electron_up, electron_down, electron_up, electron_down,
//            beryllium_nucleus};
//    const zsvm::SphericalECGContext context =
//            zsvm::SphericalECGContext::create(particles, 3);
//
//    const std::size_t basis_size = 100;
//    std::vector<std::vector<double>> basis;
//    for (std::size_t i = 0; i < basis_size; ++i) {
//        basis.push_back(context.random_correlation_coefficients());
//    }
//
//    std::vector<Eigen::MatrixXd> basis_matrices;
//    for (const auto &correlation_coefficients : basis) {
//        basis_matrices.push_back(
//                context.gaussian_parameter_matrix(correlation_coefficients));
//    }
//
//    zsvm::RealVariationalSolver solver;
//    solver.set_basis_size_destructive(basis_size);
//
//    for (std::size_t i = 0; i < basis_size; ++i) {
//        for (std::size_t j = 0; j < basis_size; ++j) {
//            double overlap_matrix_element, hamiltonian_matrix_element;
//            context.compute_matrix_elements(
//                    overlap_matrix_element, hamiltonian_matrix_element,
//                    basis_matrices[i], basis_matrices[j]);
//            solver.set_overlap_matrix_element(
//                    i, j, overlap_matrix_element);
//            solver.set_hamiltonian_matrix_element(
//                    i, j, hamiltonian_matrix_element);
//        }
//    }
//
//    PRINT_EXECUTION_TIME(std::cout << solver.get_eigenvalue(0) << std::endl);

    return 0;
}
