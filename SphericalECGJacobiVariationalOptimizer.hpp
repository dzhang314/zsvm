#ifndef ZSVM_SPHERICAL_ECG_JACOBI_VARIATIONAL_OPTIMIZER_HPP_INCLUDED
#define ZSVM_SPHERICAL_ECG_JACOBI_VARIATIONAL_OPTIMIZER_HPP_INCLUDED

// C++ standard library headers
#include <chrono>
#include <cstddef> // for std::size_t
#include <random> // for std::random_device
#include <string>
#include <vector>

// OpenMP multithreading headers
#include <omp.h>

// Boost library headers
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/random.hpp> // for boost::random::independent_bits_engine etc.

// Project-specific headers
#include "AmoebaOptimizer.hpp"
#include "RealVariationalSolver.hpp"
#include "SphericalECGJacobiContext.hpp"

namespace zsvm {

    template <typename T>
    class SphericalECGJacobiVariationalOptimizer {

        typedef boost::random::mt19937 random_generator_t;

        typedef std::array<
                random_generator_t::result_type,
                random_generator_t::state_size> random_state_t;

        typedef boost::random::independent_bits_engine<
                random_generator_t,
                std::numeric_limits<T>::digits,
                boost::multiprecision::cpp_int> random_engine_t;

    private: // ================================================================

        const std::size_t num_particles;
        const std::size_t num_parameters;
        SphericalECGJacobiContext<T> context;
        RealVariationalSolver<T> solver;
        std::vector<std::vector<T>> basis;
        std::vector<std::vector<T>> basis_matrices;
        random_engine_t random_engine;
        boost::random::normal_distribution<T> correlation_distribution;

        static random_engine_t properly_seeded_random_engine() {
            random_state_t seed_data;
            std::random_device nondet_random_source;
            std::generate(seed_data.begin(), seed_data.end(),
                          std::ref(nondet_random_source));
            seed_data[0] = static_cast<random_generator_t::result_type>(
                    std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::system_clock::now().time_since_epoch()
                    ).count() % std::numeric_limits<
                            random_generator_t::result_type>::max());
            boost::random::seed_seq seed_sequence(seed_data.begin(),
                                                  seed_data.end());
            random_engine_t random_engine(seed_sequence);
            return random_engine;
        }

    public: // =================================================================

        explicit SphericalECGJacobiVariationalOptimizer(
                long long int space_dimension,
                const std::vector<zsvm::Particle<T>> &particles,
                const std::string &mass_carrier,
                const std::string &charge_carrier)
                : num_particles(particles.size()),
                  num_parameters(num_particles * (num_particles - 1) / 2),
                  context(SphericalECGJacobiContext<T>::create(
                          particles, mass_carrier, charge_carrier,
                          space_dimension)),
                  random_engine(properly_seeded_random_engine()),
                  correlation_distribution(0, 3) {}

        T get_ground_state_energy() {
            return solver.get_eigenvalue(0);
        }

    public: // =================================================================

        const std::vector<std::vector<T>> &get_basis() {
            return basis;
        }

        void set_basis(const std::vector<std::vector<T>> &input_basis) {
            basis = input_basis;
            recompute_basis_matrices();
            recompute_solver_matrices();
        }

        T augmented_ground_state_energy(
                const T *__restrict__ new_basis_matrix) {
            if (solver.empty()) {
                T overlap_element, hamiltonian_element;
                context.evaluate_matrix_elements(
                        overlap_element, hamiltonian_element,
                        new_basis_matrix, new_basis_matrix);
                return hamiltonian_element / overlap_element;
            }
            const std::size_t basis_size = basis_matrices.size();
            std::vector<T> new_overlap_column(basis_size + 1);
            std::vector<T> new_hamiltonian_column(basis_size + 1);
            for (std::size_t i = 0; i < basis_size; ++i) {
                context.evaluate_matrix_elements(
                        new_overlap_column[i], new_hamiltonian_column[i],
                        basis_matrices[i].data(), new_basis_matrix);
            }
            context.evaluate_matrix_elements(
                    new_overlap_column[basis_size],
                    new_hamiltonian_column[basis_size],
                    new_basis_matrix, new_basis_matrix);
            return solver.minimum_augmented_eigenvalue(
                    new_overlap_column.data(),
                    new_hamiltonian_column.data());
        }

        void random_basis_element(T *__restrict__ basis_element) {
            for (std::size_t i = 0; i < num_parameters; ++i) {
                basis_element[i] = correlation_distribution(random_engine);
            }
        }

        void construct_basis_matrix(
                T *__restrict__ basis_matrix,
                const T *__restrict__ basis_element) {
            context.gaussian_parameter_matrix(basis_element, basis_matrix);
        }

        void add_basis_element(const std::vector<T> &basis_element) {
            basis.push_back(basis_element);
            std::vector<T> basis_matrix(num_parameters);
            construct_basis_matrix(basis_matrix.data(), basis_element.data());
            basis_matrices.push_back(basis_matrix);
            const std::size_t basis_size = basis.size();
            solver.set_basis_size_conservative(basis_size);
            for (std::size_t i = 0; i < basis_size - 1; ++i) {
                context.evaluate_matrix_elements(
                        solver.overlap_matrix_element(i, basis_size - 1),
                        solver.hamiltonian_matrix_element(i, basis_size - 1),
                        basis_matrices[i].data(),
                        basis_matrices[basis_size - 1].data());
                solver.overlap_matrix_element(basis_size - 1, i) =
                        solver.overlap_matrix_element(i, basis_size - 1);
                solver.hamiltonian_matrix_element(basis_size - 1, i) =
                        solver.hamiltonian_matrix_element(i, basis_size - 1);
            }
            context.evaluate_matrix_elements(
                    solver.overlap_matrix_element(
                            basis_size - 1, basis_size - 1),
                    solver.hamiltonian_matrix_element(
                            basis_size - 1, basis_size - 1),
                    basis_matrices[basis_size - 1].data(),
                    basis_matrices[basis_size - 1].data());
        }

        void replace_basis_element(const std::vector<T> &basis_element) {
            const std::size_t basis_size = basis.size();
            basis[basis_size - 1] = basis_element;
            construct_basis_matrix(basis_matrices[basis_size - 1].data(),
                                   basis_element.data());
            for (std::size_t i = 0; i < basis_size - 1; ++i) {
                context.evaluate_matrix_elements(
                        solver.overlap_matrix_element(i, basis_size - 1),
                        solver.hamiltonian_matrix_element(i, basis_size - 1),
                        basis_matrices[i].data(),
                        basis_matrices[basis_size - 1].data());
                solver.overlap_matrix_element(basis_size - 1, i) =
                        solver.overlap_matrix_element(i, basis_size - 1);
                solver.hamiltonian_matrix_element(basis_size - 1, i) =
                        solver.hamiltonian_matrix_element(i, basis_size - 1);
            }
            context.evaluate_matrix_elements(
                    solver.overlap_matrix_element(
                            basis_size - 1, basis_size - 1),
                    solver.hamiltonian_matrix_element(
                            basis_size - 1, basis_size - 1),
                    basis_matrices[basis_size - 1].data(),
                    basis_matrices[basis_size - 1].data());
        }

        std::vector<T> last_basis_element() {
            return basis.back();
        }

        bool expand_amoeba(std::size_t num_trials,
                           const T &initial_step_size,
                           std::size_t max_steps) {
            std::vector<T> new_basis_element(num_parameters);
            std::vector<T> best_basis_element;
            T best_energy = solver.empty()
                            ? std::numeric_limits<T>::max()
                            : solver.get_eigenvalue(0);
            for (std::size_t trial = 0; trial < num_trials; ++trial) {
                random_basis_element(new_basis_element.data());
                const T new_energy = refine_amoeba(
                        new_basis_element.data(), initial_step_size, max_steps);
                if (new_energy < best_energy) {
                    best_basis_element = new_basis_element;
                    best_energy = new_energy;
                }
            }
            if (!best_basis_element.empty()) {
                add_basis_element(best_basis_element);
                return true;
            } else {
                return false;
            }
        }

    public: // =================================================================

        T refine_amoeba(T *__restrict__ basis_element,
                        const T &initial_step_size,
                        std::size_t max_steps) {
            std::vector<T> basis_matrix(num_parameters);
            dznl::AmoebaOptimizer<T> amoeba(
                    basis_element, num_parameters, initial_step_size,
                    [&](const T *b) {
                        construct_basis_matrix(basis_matrix.data(), b);
                        return augmented_ground_state_energy(
                                basis_matrix.data());
                    });
            for (std::size_t step = 0; step < max_steps; ++step) {
                amoeba.step();
            }
            return amoeba.current_minimum(basis_element);
        }

        void recompute_basis_matrices() {
            basis_matrices.clear();
            for (const auto &basis_element : basis) {
                std::vector<T> basis_matrix(num_parameters);
                construct_basis_matrix(basis_matrix.data(),
                                       basis_element.data());
                basis_matrices.push_back(basis_matrix);
            }
        }

        void recompute_solver_matrices() {
            const std::size_t num_basis_matrices = basis_matrices.size();
            solver.set_basis_size_destructive(num_basis_matrices);
            for (std::size_t i = 0; i < num_basis_matrices; ++i) {
                for (std::size_t j = 0; j < num_basis_matrices; ++j) {
                    context.evaluate_matrix_elements(
                            solver.overlap_matrix_element(i, j),
                            solver.hamiltonian_matrix_element(i, j),
                            basis_matrices[i].data(),
                            basis_matrices[j].data());
                }
            }
        }

    }; // class SphericalECGJacobiVariationalOptimizer

} // namespace zsvm

#endif // ZSVM_SPHERICAL_ECG_JACOBI_VARIATIONAL_OPTIMIZER_HPP_INCLUDED
