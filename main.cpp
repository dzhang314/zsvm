// C++ standard library headers
#include <iostream>
#include <vector>

// Project-specific headers
#include "Particle.hpp"
#include "RealSolver.hpp"

int main() {
    const zsvm::Particle electron_up = {0, 1.0, -1.0, zsvm::Spin::UP};
    const zsvm::Particle electron_down = {0, 1.0, -1.0, zsvm::Spin::DOWN};
    const zsvm::Particle positron_up = {1, 1.0, +1.0, zsvm::Spin::UP};
    const zsvm::Particle positron_down = {1, 1.0, +1.0, zsvm::Spin::DOWN};
    std::vector<zsvm::Particle> particles = {
            electron_up, electron_down, positron_up, positron_down};
    zsvm::RealSolver solver(particles);
    return 0;
}
