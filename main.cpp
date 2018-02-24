// C++ standard library headers
#include <chrono>
#include <iostream>
#include <vector>

// Eigen linear algebra library headers
#include <Eigen/Eigenvalues>

// Project-specific headers
#include "Particle.hpp"
#include "SphericalECGContext.hpp"


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

    std::vector<double> v = context.random_correlation_coefficients();
    std::vector<double> w = context.random_correlation_coefficients();
    Eigen::MatrixXd a = context.gaussian_parameter_matrix(v);
    Eigen::MatrixXd b = context.gaussian_parameter_matrix(w);

    {
        double overlap_matrix_element, hamiltonian_matrix_element;
        std::chrono::high_resolution_clock::time_point start =
                std::chrono::high_resolution_clock::now();
        context.matrix_elements(
                overlap_matrix_element, hamiltonian_matrix_element, a, b);
        std::chrono::high_resolution_clock::time_point stop =
                std::chrono::high_resolution_clock::now();
        std::cout << overlap_matrix_element << std::endl;
        std::cout << hamiltonian_matrix_element << std::endl;
        std::cout << std::chrono::duration_cast<std::chrono::microseconds>(
                stop - start).count() << std::endl;
    }

//    const std::size_t basis_size = 100;
//    std::vector<std::vector<double>> basis;
//    for (std::size_t i = 0; i < basis_size; ++i) {
//        basis.push_back(solver.random_correlation_coefficients());
//    }
//
//    std::vector<Eigen::MatrixXd> basis_matrices;
//    for (const auto &correlation_coefficients : basis) {
//        basis_matrices.push_back(
//                solver.gaussian_parameter_matrix(correlation_coefficients));
//    }
//
//    Eigen::MatrixXd hamiltonian_matrix(basis_size, basis_size);
//    Eigen::MatrixXd overlap_matrix(basis_size, basis_size);
//
//    for (std::size_t i = 0; i < basis_size; ++i) {
//        for (std::size_t j = 0; j < basis_size; ++j) {
//            hamiltonian_matrix(i, j) =
//                    solver.kinetic_matrix_element(
//                            basis_matrices[i], basis_matrices[j]) +
//                    solver.coulomb_matrix_element(
//                            basis_matrices[i], basis_matrices[j]);
//            overlap_matrix(i, j) =
//                    solver.overlap_matrix_element(
//                            basis_matrices[i], basis_matrices[j]);
//        }
//    }
//
//    Eigen::GeneralizedSelfAdjointEigenSolver energy_solver(
//            hamiltonian_matrix, overlap_matrix,
//            Eigen::EigenvaluesOnly | Eigen::Ax_lBx);
//    double best_energy = energy_solver.eigenvalues()[0];

//    while (true) {
//        for (std::size_t i = 0; i < basis_size; ++i) {
//            std::vector<double> old_basis_element = basis[i];
//            basis[i] = solver.random_correlation_coefficients();
//            basis_matrices[i] = solver.gaussian_parameter_matrix(basis[i]);
//            for (std::size_t j = 0; j < basis_size; ++j) {
//                hamiltonian_matrix(i, j) =
//                        solver.kinetic_matrix_element(
//                                basis_matrices[i], basis_matrices[j]) +
//                        solver.coulomb_matrix_element(
//                                basis_matrices[i], basis_matrices[j]);
//                hamiltonian_matrix(j, i) =
//                        solver.kinetic_matrix_element(
//                                basis_matrices[j], basis_matrices[i]) +
//                        solver.coulomb_matrix_element(
//                                basis_matrices[j], basis_matrices[i]);
//                overlap_matrix(i, j) =
//                        solver.overlap_matrix_element(
//                                basis_matrices[i], basis_matrices[j]);
//                overlap_matrix(j, i) =
//                        solver.overlap_matrix_element(
//                                basis_matrices[j], basis_matrices[i]);
//            }
//            Eigen::GeneralizedSelfAdjointEigenSolver test_energy_solver(
//                    hamiltonian_matrix, overlap_matrix,
//                    Eigen::EigenvaluesOnly | Eigen::Ax_lBx);
//            double test_energy = test_energy_solver.eigenvalues()[0];
//            if (test_energy < best_energy) {
//                best_energy = test_energy;
//                std::cout << best_energy << std::endl;
//            } else {
//                basis[i] = old_basis_element;
//                basis_matrices[i] = solver.gaussian_parameter_matrix(basis[i]);
//                for (std::size_t j = 0; j < basis_size; ++j) {
//                    hamiltonian_matrix(i, j) =
//                            solver.kinetic_matrix_element(
//                                    basis_matrices[i], basis_matrices[j]) +
//                            solver.coulomb_matrix_element(
//                                    basis_matrices[i], basis_matrices[j]);
//                    hamiltonian_matrix(j, i) =
//                            solver.kinetic_matrix_element(
//                                    basis_matrices[j], basis_matrices[i]) +
//                            solver.coulomb_matrix_element(
//                                    basis_matrices[j], basis_matrices[i]);
//                    overlap_matrix(i, j) =
//                            solver.overlap_matrix_element(
//                                    basis_matrices[i], basis_matrices[j]);
//                    overlap_matrix(j, i) =
//                            solver.overlap_matrix_element(
//                                    basis_matrices[j], basis_matrices[i]);
//                }
//            }
//        }
//    }

    return 0;
}
