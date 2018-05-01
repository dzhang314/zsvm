// C++ standard library headers
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
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

typedef boost::multiprecision::number<
        boost::multiprecision::cpp_bin_float<30>> bigfloat_t;

//typedef double bigfloat_t;

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
    particles.emplace_back(0, -1, electron_carriers);
    particles.emplace_back(1, +1, positron_carriers);
    // particles.emplace_back(1, -1, positron_carriers);

    std::size_t num_threads = 1;
#pragma omp parallel
    {
        if (omp_get_thread_num() == 0) {
            num_threads = static_cast<std::size_t>(omp_get_num_threads());
        }
    }

    if (num_threads == 1) {
        std::cout << "Performing serial basis expansion.\n";
    } else {
        std::cout << "Performing basis expansion using "
                  << num_threads << " threads.\n";
    }

    const std::size_t target_size = 100'000;
    const std::size_t trials_per_thread = 1;
    const std::size_t max_iterations = 10;

    std::ifstream basis_input_file("basis_0000000600.txt");
    std::vector<std::string> basis_input_lines;
    for (std::string line; std::getline(basis_input_file, line);) {
        basis_input_lines.push_back(line);
    }

    std::vector<std::vector<std::string>> basis_input_items;
    for (const auto &line : basis_input_lines) {
        std::stringstream line_stream(line);
        basis_input_items.emplace_back();
        for (std::string item; std::getline(line_stream, item, '\t');) {
            if (!item.empty()) {
                basis_input_items.back().push_back(item);
            }
        }
    }

    basis_input_items.erase(
            std::remove_if(basis_input_items.begin(), basis_input_items.end(),
                           [](const std::vector<std::string> &items) {
                               return items.empty();
                           }),
            basis_input_items.end());

//    for (const auto &line : basis_input_items) {
//        std::cout << "LINE(";
//        for (const auto &item : line) {
//            std::cout << "ITEM(" << item << "), ";
//        }
//        std::cout << ")\n";
//    }
//    std::cout << std::endl;

    std::vector<std::vector<bigfloat_t>> input_basis;
    for (const auto &line : basis_input_items) {
        input_basis.emplace_back();
        for (const auto &item : line) {
            std::stringstream item_stream(item);
            bigfloat_t x;
            item_stream >> x;
            input_basis.back().push_back(x);
        }
    }

//    for (const auto &basis_element : input_basis) {
//        std::cout << "LINE(";
//        for (const auto &coefficient : basis_element) {
//            std::cout << "NUMBER(" << coefficient << "), ";
//        }
//        std::cout << ")\n";
//    }
//    std::cout << std::endl;

    std::chrono::high_resolution_clock::time_point start =
            std::chrono::high_resolution_clock::now();

    std::vector<bigfloat_t> best_basis_element;
    std::vector<std::vector<bigfloat_t>> thread_basis_elements(num_threads);
    std::vector<bigfloat_t> thread_energies(num_threads);

    std::size_t min_index;
    bool any_found;

#pragma omp parallel
    {
        const int thread_num = omp_get_thread_num();
        zsvm::SphericalECGJacobiVariationalOptimizer<bigfloat_t> optimizer(
                3, particles, mass_carrier, charge_carrier);
        optimizer.set_basis(input_basis);

        for (std::size_t basis_size = optimizer.get_basis().size();
             basis_size <= target_size;
             basis_size = optimizer.get_basis().size()) {
            bool found = optimizer.expand_amoeba(
                    trials_per_thread * num_threads, 1, max_iterations);
            if (found) {
                thread_basis_elements[thread_num] =
                        optimizer.last_basis_element();
                thread_energies[thread_num] =
                        optimizer.get_ground_state_energy();
            } else {
                thread_energies[thread_num] =
                        std::numeric_limits<bigfloat_t>::max();
            }
#pragma omp barrier
            if (thread_num == 0) {
                min_index = static_cast<std::size_t>(
                        std::min_element(thread_energies.begin(),
                                         thread_energies.end()) -
                        thread_energies.begin());
                if (thread_energies[min_index] ==
                    std::numeric_limits<bigfloat_t>::max()) {
                    any_found = false;
                    std::cout << basis_size << '\t'
                              << "FAILED TO IMPROVE ENERGY" << std::endl;
                } else {
                    any_found = true;
                    std::cout << basis_size << '\t'
                              << thread_energies[min_index] << std::endl;
                }
            }
#pragma omp barrier
            if (found && any_found) {
                optimizer.replace_basis_element(
                        thread_basis_elements[min_index]);
            } else if (!found && any_found) {
                optimizer.add_basis_element(
                        thread_basis_elements[min_index]);
            }
#pragma omp barrier
            if (thread_num == 0 && any_found) {
                std::ostringstream file_name;
                file_name << "basis_" << std::setfill('0') << std::setw(10)
                          << basis_size << ".txt";
                std::ofstream basis_file(file_name.str());
                basis_file << std::scientific << std::showpos
                           << std::setprecision(std::numeric_limits<
                                   bigfloat_t>::max_digits10 - 1);
                for (const auto &basis_element : optimizer.get_basis()) {
                    for (std::size_t i = 0; i < basis_element.size(); ++i) {
                        if (i > 0) { basis_file << '\t'; }
                        basis_file << basis_element[i];
                    }
                    basis_file << '\n';
                }
            }
        }
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
