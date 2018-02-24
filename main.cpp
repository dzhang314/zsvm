// C++ standard library headers
#include <iostream>
#include <vector>

// Project-specific headers
#include "Particle.hpp"
#include "RealSolver.hpp"

int main() {
    const zsvm::Particle electron_up = {0, 1.0, -1.0, zsvm::Spin::UP};
    const zsvm::Particle electron_down = {0, 1.0, -1.0, zsvm::Spin::DOWN};
    const zsvm::Particle beryllium_nucleus = {1, 16538.028978017737,
                                              +4.0, zsvm::Spin::UP};
    std::vector<zsvm::Particle> particles = {
            electron_up, electron_down, electron_up, electron_down,
            beryllium_nucleus};
    zsvm::RealSolver solver(particles, 3);
    std::vector<double> v = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 1.0};
    std::vector<double> w = {1.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0};
    Eigen::MatrixXd a = solver.gaussian_parameter_matrix(v);
    Eigen::MatrixXd b = solver.gaussian_parameter_matrix(w);
    std::cout << solver.overlap_matrix_element(a, b) << std::endl;
    std::cout << solver.overlap_matrix_element(b, a) << std::endl;
    std::cout << solver.overlap_matrix_element(a, b) -
                 solver.overlap_matrix_element(b, a) << std::endl;
    std::cout << solver.kinetic_matrix_element(a, b) << std::endl;
    std::cout << solver.kinetic_matrix_element(b, a) << std::endl;
    std::cout << solver.kinetic_matrix_element(a, b) -
                 solver.kinetic_matrix_element(b, a) << std::endl;
    std::cout << solver.coulomb_matrix_element(a, b) << std::endl;
    std::cout << solver.coulomb_matrix_element(b, a) << std::endl;
    std::cout << solver.coulomb_matrix_element(a, b) -
                 solver.coulomb_matrix_element(b, a) << std::endl;
    return 0;
}
