// C++ standard library headers
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

// OpenMP multithreading headers
#include <omp.h>

// Boost library headers
#include <boost/io/ios_state.hpp>

#define BOOST_MP_MIN_EXPONENT_DIGITS 3

#include <boost/multiprecision/cpp_bin_float.hpp>

// Project-specific headers
#include "Particle.hpp"
#include "SphericalECGJacobiVariationalOptimizer.hpp"

typedef boost::multiprecision::cpp_bin_float_50 bigfloat_t;

int main(int, char **) {
    std::cout << std::setprecision(
            std::numeric_limits<bigfloat_t>::max_digits10 - 1);

    std::vector<zsvm::Particle<bigfloat_t>> particles;
    const std::string mass_carrier("mass");
    const std::string charge_carrier("charge");

    std::map<std::string, bigfloat_t> electron_carriers;
    electron_carriers.insert({mass_carrier, 1});
    electron_carriers.insert({charge_carrier, -1});

    std::map<std::string, bigfloat_t> positron_carriers;
    positron_carriers.insert({mass_carrier, 1});
    positron_carriers.insert({charge_carrier, +1});

    particles.emplace_back(0, +1, electron_carriers);
    // particles.emplace_back(0, -1, electron_carriers);
    particles.emplace_back(1, +1, positron_carriers);
    // particles.emplace_back(1, -1, positron_carriers);

    std::size_t num_threads = 1;
#pragma omp parallel
    {
        if (omp_get_thread_num() == 0) {
            num_threads = (std::size_t) omp_get_num_threads();
        }
    }

    if (num_threads == 1) {
        std::cout << "Performing serial basis expansion.\n";
    } else {
        std::cout << "Performing basis expansion using "
                  << num_threads << " threads.\n";
    }

    const std::size_t target_size = 2000;
    const std::size_t trials = 32;
    const std::size_t max_iterations = 50;

    std::chrono::high_resolution_clock::time_point start =
            std::chrono::high_resolution_clock::now();

    std::vector<std::vector<bigfloat_t>> best_basis;
    std::vector<std::vector<std::vector<bigfloat_t>>> thread_bases(num_threads);
    std::vector<bigfloat_t> thread_energies(num_threads);

    for (std::size_t basis_size = 1; basis_size <= target_size; ++basis_size) {
#pragma omp parallel
        {
            zsvm::SphericalECGJacobiVariationalOptimizer<bigfloat_t> optimizer(
                    3, particles, mass_carrier, charge_carrier);
            optimizer.basis = best_basis;
            optimizer.recompute_basis_matrices();
            optimizer.recompute_solver_matrices();
            optimizer.expand_amoeba(trials / num_threads, 1.0, max_iterations);
            thread_bases[omp_get_thread_num()] = optimizer.basis;
            thread_energies[omp_get_thread_num()] =
                    optimizer.get_ground_state_energy();
        }
        std::size_t min_index = static_cast<std::size_t>(
                std::min_element(thread_energies.begin(),
                                 thread_energies.end()) -
                thread_energies.begin());
        best_basis = thread_bases[min_index];
        {
            std::ostringstream file_name;
            file_name << "basis_" << std::setfill('0') << std::setw(10)
                      << basis_size << ".txt";
            std::ofstream basis_file(file_name.str());
            basis_file << std::scientific << std::showpos << std::setprecision(
                    std::numeric_limits<bigfloat_t>::max_digits10 - 1);
            for (const auto &basis_element : best_basis) {
                for (std::size_t i = 0; i < basis_element.size(); ++i) {
                    if (i > 0) { basis_file << '\t'; }
                    basis_file << basis_element[i];
                }
                basis_file << '\n';
            }
        }
        std::cout << basis_size << '\t' << thread_energies[min_index]
                  << std::endl;
    }

    std::chrono::high_resolution_clock::time_point stop =
            std::chrono::high_resolution_clock::now();
    std::cout << "Time elapsed (seconds): "
              << std::chrono::duration_cast<
                      std::chrono::milliseconds>(
                      stop - start).count() / 1000.0 << std::endl;

//    if (argc != 2) {
//        std::cout << "Usage: " << argv[0] << " SCRIPT_FILE" << std::endl;
//        return 1;
//    }
//    std::string script_file_name(argv[1]);
//    {
//        std::ifstream script_file(script_file_name);
//        if (!script_file.good()) {
//            std::cout << "Error: could not open script file "
//                      << script_file_name << "." << std::endl;
//            return 2;
//        }
//    }
//    std::cout << std::scientific;
//    zsvm::ScriptInterpreter interpreter(script_file_name);
//    interpreter.run();
//    return 0;
}
