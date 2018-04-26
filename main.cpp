// C++ standard library headers
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>
#include <random>

#include <omp.h>

// Project-specific headers
#include "EnumTypes.hpp"
#include "Particle.hpp"
#include "RealVariationalSolver.hpp"
#include "AmoebaOptimizer.hpp"
#include "ScriptParser.hpp"
#include "SphericalECGOverlapContext.hpp"
#include "SphericalECGJacobiContext.hpp"

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

    public:

        const std::size_t num_particles;
        const std::size_t num_parameters;
        SphericalECGOverlapContext context;
        RealVariationalSolver solver;
        std::vector<std::vector<double>> basis;
        std::vector<std::vector<double>> basis_matrices;

        explicit SphericalECGVariationalOptimizer(
                long long int space_dimension,
                const std::vector<Particle<double>> &particles,
                const std::map<std::string, DispersionRelation<double>> &
                dispersion_relations,
                const std::map<std::string, ConfiningPotential<double>> &
                confining_potentials,
                const std::map<std::string, PairwisePotential<double>> &
                pairwise_potentials)
                : num_particles(particles.size()),
                  num_parameters(num_particles * (num_particles + 1) / 2),
                  context(space_dimension, particles, dispersion_relations,
                          confining_potentials, pairwise_potentials) {}

        double get_ground_state_energy() {
            return solver.get_eigenvalue(0);
        }

    public: // =================================================================

        double augmented_ground_state_energy(
                const double *__restrict__ new_basis_matrix) {
            if (solver.empty()) {
                double overlap_element, hamiltonian_element;
                context.evaluate_matrix_elements(
                        overlap_element, hamiltonian_element,
                        new_basis_matrix, new_basis_matrix);
                return hamiltonian_element / overlap_element;
            } else {
                const std::size_t basis_size = basis_matrices.size();
                std::vector<double> new_overlap_column(basis_size + 1);
                std::vector<double> new_hamiltonian_column(basis_size + 1);
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
        }

        static std::mt19937_64 properly_seeded_random_engine() {
            std::array<std::mt19937_64::result_type,
                    std::mt19937_64::state_size> random_data;
            std::random_device source;
            std::generate(random_data.begin(), random_data.end(),
                          std::ref(source));
            random_data[0] = static_cast<std::mt19937_64::result_type>(
                    std::chrono::duration_cast<std::chrono::milliseconds>(
                            std::chrono::system_clock::now().time_since_epoch()
                    ).count());
            std::seed_seq seeds(random_data.begin(), random_data.end());
            return std::mt19937_64(seeds);
        }

        void random_basis_element(double *__restrict__ basis_element) {
            static std::mt19937_64 random_engine = properly_seeded_random_engine();
            static std::normal_distribution<double>
                    correlation_distribution(0.0, 3.0);
            for (std::size_t i = 0; i < num_parameters; ++i) {
                basis_element[i] = correlation_distribution(random_engine);
            }
        }

        void construct_basis_matrix(
                double *__restrict__ basis_matrix,
                const double *__restrict__ basis_element) {
            for (std::size_t i = 0, k = 0; i < num_particles; ++i) {
                for (std::size_t j = i; j < num_particles; ++j, ++k) {
                    const double x = std::exp(basis_element[k]);
                    basis_matrix[i * (i + 3) / 2] += x;
                    basis_matrix[j * (j + 3) / 2] += x;
                    basis_matrix[j * (j + 1) / 2 + i] -= x;
                }
            }
        }

        void expand_random(std::size_t num_trials) {
            std::vector<double> new_basis_element(num_parameters);
            std::vector<double> new_basis_matrix(num_parameters);
            std::vector<double> best_basis_element;
            std::vector<double> best_basis_matrix;
            double best_energy = solver.empty()
                                 ? std::numeric_limits<double>::max()
                                 : solver.get_eigenvalue(0);
            for (std::size_t trial = 0; trial < num_trials; ++trial) {
                random_basis_element(new_basis_element.data());
                construct_basis_matrix(new_basis_matrix.data(),
                                       new_basis_element.data());
                const double new_energy =
                        augmented_ground_state_energy(new_basis_matrix.data());
                if (new_energy < best_energy) {
                    best_basis_element = new_basis_element;
                    best_basis_matrix = new_basis_matrix;
                    best_energy = new_energy;
                }
            }
            if (!best_basis_element.empty() && !best_basis_matrix.empty()) {
                basis.push_back(best_basis_element);
                basis_matrices.push_back(best_basis_matrix);
                recompute_solver_matrices();
            }
        }

        void expand_amoeba(std::size_t num_trials,
                           double initial_step_size,
                           std::size_t max_steps) {
            std::vector<double> new_basis_element(num_parameters);
            std::vector<double> best_basis_element;
            double best_energy = solver.empty()
                                 ? std::numeric_limits<double>::max()
                                 : solver.get_eigenvalue(0);
            for (std::size_t trial = 0; trial < num_trials; ++trial) {
                random_basis_element(new_basis_element.data());
                const double new_energy = refine_amoeba(
                        new_basis_element.data(), initial_step_size, max_steps);
                if (new_energy < best_energy) {
                    best_basis_element = new_basis_element;
                    best_energy = new_energy;
                }
            }
            if (!best_basis_element.empty()) {
                basis.push_back(best_basis_element);
                std::vector<double> best_basis_matrix(num_parameters);
                construct_basis_matrix(best_basis_matrix.data(),
                                       best_basis_element.data());
                basis_matrices.push_back(best_basis_matrix);
                recompute_solver_matrices();
            }
        }

    public: // =================================================================

        double refine_amoeba(double *__restrict__ basis_element,
                             double initial_step_size,
                             std::size_t max_steps) {
            std::vector<double> basis_matrix(num_parameters);
            dznl::AmoebaOptimizer amoeba(
                    basis_element, num_parameters, initial_step_size,
                    [&](const double *b) {
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
                std::vector<double> basis_matrix(num_parameters);
                construct_basis_matrix(basis_matrix.data(),
                                       basis_element.data());
                basis_matrices.push_back(basis_matrix);
            }
        }

        void recompute_solver_matrices() {
            solver.set_basis_size_destructive(basis_matrices.size());
            for (std::size_t i = 0; i < basis_matrices.size(); ++i) {
                for (std::size_t j = 0; j < basis_matrices.size(); ++j) {
                    context.evaluate_matrix_elements(
                            solver.overlap_matrix_element(i, j),
                            solver.hamiltonian_matrix_element(i, j),
                            basis_matrices[i].data(),
                            basis_matrices[j].data());
                }
            }
        }

    }; // class SphericalECGVariationalOptimizer

} // namespace zsvm


namespace zsvm {

    class SphericalECGJacobiVariationalOptimizer {

    public:

        const std::size_t num_particles;
        const std::size_t num_parameters;
        SphericalECGJacobiContext<double> context;
        RealVariationalSolver solver;
        std::vector<std::vector<double>> basis;
        std::vector<std::vector<double>> basis_matrices;

        explicit SphericalECGJacobiVariationalOptimizer(
                long long int space_dimension,
                const std::vector<Particle<double>> &particles,
                const std::string &mass_carrier,
                const std::string &charge_carrier)
                : num_particles(particles.size()),
                  num_parameters(num_particles * (num_particles - 1) / 2),
                  context(SphericalECGJacobiContext<double>::create(
                          particles, mass_carrier, charge_carrier,
                          space_dimension)) {}

        double get_ground_state_energy() {
            return solver.get_eigenvalue(0);
        }

    public: // =================================================================

        double augmented_ground_state_energy(
                const double *__restrict__ new_basis_matrix) {
            if (solver.empty()) {
                double overlap_element, hamiltonian_element;
                context.evaluate_matrix_elements(
                        overlap_element, hamiltonian_element,
                        new_basis_matrix, new_basis_matrix);
                return hamiltonian_element / overlap_element;
            } else {
                const std::size_t basis_size = basis_matrices.size();
                std::vector<double> new_overlap_column(basis_size + 1);
                std::vector<double> new_hamiltonian_column(basis_size + 1);
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
        }

        static std::mt19937_64 properly_seeded_random_engine() {
            std::array<std::mt19937_64::result_type,
                    std::mt19937_64::state_size> random_data;
            std::random_device source;
            std::generate(random_data.begin(), random_data.end(),
                          std::ref(source));
            random_data[0] = static_cast<std::mt19937_64::result_type>(
                    std::chrono::duration_cast<std::chrono::milliseconds>(
                            std::chrono::system_clock::now().time_since_epoch()
                    ).count());
            std::seed_seq seeds(random_data.begin(), random_data.end());
            return std::mt19937_64(seeds);
        }

        void random_basis_element(double *__restrict__ basis_element) {
            static std::mt19937_64 random_engine = properly_seeded_random_engine();
            static std::normal_distribution<double>
                    correlation_distribution(0.0, 3.0);
            for (std::size_t i = 0; i < num_parameters; ++i) {
                basis_element[i] = correlation_distribution(random_engine);
            }
        }

        void construct_basis_matrix(
                double *__restrict__ basis_matrix,
                const double *__restrict__ basis_element) {
            context.gaussian_parameter_matrix(basis_element, basis_matrix);
        }

        void expand_amoeba(std::size_t num_trials,
                           double initial_step_size,
                           std::size_t max_steps) {
            std::vector<double> new_basis_element(num_parameters);
            std::vector<double> best_basis_element;
            double best_energy = solver.empty()
                                 ? std::numeric_limits<double>::max()
                                 : solver.get_eigenvalue(0);
            for (std::size_t trial = 0; trial < num_trials; ++trial) {
                random_basis_element(new_basis_element.data());
                const double new_energy = refine_amoeba(
                        new_basis_element.data(), initial_step_size, max_steps);
                if (new_energy < best_energy) {
                    best_basis_element = new_basis_element;
                    best_energy = new_energy;
                }
            }
            if (!best_basis_element.empty()) {
                basis.push_back(best_basis_element);
                std::vector<double> best_basis_matrix(num_parameters);
                construct_basis_matrix(best_basis_matrix.data(),
                                       best_basis_element.data());
                basis_matrices.push_back(best_basis_matrix);
                recompute_solver_matrices();
            }
        }

    public: // =================================================================

        double refine_amoeba(double *__restrict__ basis_element,
                             double initial_step_size,
                             std::size_t max_steps) {
            std::vector<double> basis_matrix(num_parameters);
            dznl::AmoebaOptimizer amoeba(
                    basis_element, num_parameters, initial_step_size,
                    [&](const double *b) {
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
                std::vector<double> basis_matrix(num_parameters);
                construct_basis_matrix(basis_matrix.data(),
                                       basis_element.data());
                basis_matrices.push_back(basis_matrix);
            }
        }

        void recompute_solver_matrices() {
            solver.set_basis_size_destructive(basis_matrices.size());
            for (std::size_t i = 0; i < basis_matrices.size(); ++i) {
                for (std::size_t j = 0; j < basis_matrices.size(); ++j) {
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


namespace zsvm {

    void error_exit(const char *message) {
        std::cerr << message << std::endl;
        std::exit(EXIT_FAILURE);
    }


    void error_exit(const char *message, const std::string &arg) {
        std::fprintf(stderr, message, arg.data());
        std::exit(EXIT_FAILURE);
    }


    long long int get_integer(const ScriptCommand &command) {
        const auto &words = command.get_words();
        if (words.size() != 3
            || words[2].type != ScriptToken::Type::INTEGER) {
            std::cerr << "ERROR: Invalid command";
            // TODO: Improve reporting of this error. Don't print tokens.
            for (const auto &word : words) { std::cerr << ' ' << word; }
            std::cerr << ". Expected format is <verb> <noun> <integer>."
                      << std::endl;
            std::exit(EXIT_FAILURE);
        }
        return words[2].integer_value;
    }


    template <typename ValueType>
    ValueType get_enum_parameter(
            const ScriptCommand &command,
            const std::string &param_name,
            const std::map<std::string, ValueType> &valid_values) {
        const auto &words = command.get_words();
        const auto &params = command.get_named_parameters();
        const auto token_iter = params.find(param_name);
        if (token_iter == params.end()) {
            std::cerr << "ERROR: Command";
            for (const auto &word : words) { std::cerr << ' ' << word; }
            std::cerr << " lacks required named parameter \"" << param_name
                      << "\". Valid values are:";
            for (const auto &pair : valid_values) {
                std::cerr << " \"" << pair.first << '\"';
            }
            std::cerr << std::endl;
            std::exit(EXIT_FAILURE);
        }
        const ScriptToken &token = token_iter->second;
        if (token.type != ScriptToken::Type::IDENTIFIER) {
            std::cerr << "ERROR: Command";
            for (const auto &word : words) { std::cerr << ' ' << word; }
            std::cerr << " given invalid value " << token
                      << " for required named parameter \"" << param_name
                      << "\". Valid values are:";
            for (const auto &pair : valid_values) {
                std::cerr << " \"" << pair.first << '\"';
            }
            std::cerr << std::endl;
            std::exit(EXIT_FAILURE);
        }
        const auto value_iter = valid_values.find(token.string_value);
        if (value_iter == valid_values.end()) {
            std::cerr << "ERROR: Command";
            for (const auto &word : words) { std::cerr << ' ' << word; }
            std::cerr << " given invalid value " << token
                      << " for required named parameter \"" << param_name
                      << "\". Valid values are:";
            for (const auto &pair : valid_values) {
                std::cerr << " \"" << pair.first << '\"';
            }
            std::cerr << std::endl;
            std::exit(EXIT_FAILURE);
        }
        return value_iter->second;
    }


    double get_double_parameter(
            const ScriptCommand &command, const std::string &param_name) {
        const auto &words = command.get_words();
        const auto &params = command.get_named_parameters();
        const auto token_iter = params.find(param_name);
        if (token_iter == params.end()) {
            std::cerr << "ERROR: Command";
            for (const auto &word : words) { std::cerr << ' ' << word; }
            std::cerr << " lacks required named parameter \"" << param_name
                      << "\", which is expected to be an integer or real "
                      << "number." << std::endl;
            std::exit(EXIT_FAILURE);
        }
        const ScriptToken &token = token_iter->second;
        if (token.type != ScriptToken::Type::INTEGER
            && token.type != ScriptToken::Type::DECIMAL) {
            std::cerr << "ERROR: Command";
            for (const auto &word : words) { std::cerr << ' ' << word; }
            std::cerr << " given invalid value " << token
                      << " for required named parameter \"" << param_name
                      << "\", which is expected to be an integer or real "
                      << "number." << std::endl;
            std::exit(EXIT_FAILURE);
        }
        return token.double_value;
    }


    std::size_t get_unsigned_parameter(
            const ScriptCommand &command, const std::string &param_name) {
        const auto &words = command.get_words();
        const auto &params = command.get_named_parameters();
        const auto token_iter = params.find(param_name);
        if (token_iter == params.end()) {
            std::cerr << "ERROR: Command";
            for (const auto &word : words) { std::cerr << ' ' << word; }
            std::cerr << " lacks required named parameter \"" << param_name
                      << "\", which is expected to be an non-negative integer."
                      << std::endl;
            std::exit(EXIT_FAILURE);
        }
        const ScriptToken &token = token_iter->second;
        if (token.type != ScriptToken::Type::INTEGER
            || token.integer_value < 0) {
            std::cerr << "ERROR: Command";
            for (const auto &word : words) { std::cerr << ' ' << word; }
            std::cerr << " given invalid value " << token
                      << " for required named parameter \"" << param_name
                      << "\", which is expected to be am non-negative integer."
                      << std::endl;
            std::exit(EXIT_FAILURE);
        }
        return static_cast<std::size_t>(token.integer_value);
    }


    const std::string &get_string_parameter(
            const ScriptCommand &command, const std::string &param_name) {
        const auto &words = command.get_words();
        const auto &params = command.get_named_parameters();
        const auto token_iter = params.find(param_name);
        if (token_iter == params.end()) {
            std::cerr << "ERROR: Command";
            for (const auto &word : words) { std::cerr << ' ' << word; }
            std::cerr << " lacks required named parameter \"" << param_name
                      << "\", which is expected to be a valid identifier."
                      << std::endl;
            std::exit(EXIT_FAILURE);
        }
        const ScriptToken &token = token_iter->second;
        if (token.type != ScriptToken::Type::IDENTIFIER) {
            std::cerr << "ERROR: Command";
            for (const auto &word : words) { std::cerr << ' ' << word; }
            std::cerr << " given invalid value " << token
                      << " for required named parameter \"" << param_name
                      << "\", which is expected to be a valid identifier."
                      << std::endl;
            std::exit(EXIT_FAILURE);
        }
        return token.string_value;
    }


    const char *to_descriptive_string(ExchangeStatistics::Type type) {
        switch (type) {
            case ExchangeStatistics::Type::BOSON:
                return "bosonic exchange statistics";
            case ExchangeStatistics::Type::FERMION:
                return "fermionic exchange statistics";
            case ExchangeStatistics::Type::DISTINGUISHABLE:
                return "no exchange symmetry";
        }
        return ""; // Dummy return to suppress compiler warning.
    }


    std::string to_unsigned_spin(int spin) {
        std::ostringstream stream;
        if (spin % 2 == 0) {
            stream << spin / 2;
        } else {
            stream << spin << "/2";
        }
        return stream.str();
    }


    std::string to_signed_spin(int spin) {
        if (spin == 0) { return "0"; }
        int abs_spin = std::abs(spin);
        std::ostringstream stream;
        stream << ((spin > 0) ? '+' : '-');
        if (abs_spin % 2 == 0) {
            stream << abs_spin / 2;
        } else {
            stream << abs_spin << "/2";
        }
        return stream.str();
    }


    class ScriptInterpreter {

        ScriptParser parser;
        std::optional<long long int> space_dimension;
        std::map<std::string, std::tuple<std::size_t, ExchangeStatistics, int>>
                particle_types; // name -> (ID, exchange statistics, spin)
        std::map<std::string, DispersionRelation<double>> dispersion_relations;
        std::map<std::string, ConfiningPotential<double>> confining_potentials;
        std::map<std::string, PairwisePotential<double>> pairwise_potentials;
        std::vector<Particle<double>> particles;
        std::vector<std::vector<double>> basis;
        std::vector<std::vector<double>> basis_matrices;

    public: // ===================================================== CONSTRUCTOR

        explicit ScriptInterpreter(const std::string &script_file_name)
                : parser(script_file_name),
                  space_dimension(),
                  particle_types(),
                  dispersion_relations(),
                  confining_potentials(),
                  pairwise_potentials(),
                  particles(),
                  basis(),
                  basis_matrices() {}

    private: // ================================================================

        const std::string &get_unique_identifier(const ScriptCommand &command) {
            const auto &words = command.get_words();
            if (words.size() != 3
                || words[2].type != ScriptToken::Type::IDENTIFIER) {
                std::cerr << "ERROR: Invalid command";
                // TODO: Improve reporting of this error. Don't print tokens.
                for (const auto &word : words) { std::cerr << ' ' << word; }
                std::cerr << ". Expected format is <verb> <noun> <name> "
                          << "where <name> is a valid identifier." << std::endl;
                std::exit(EXIT_FAILURE);
            }
            const std::string &name = words[2].string_value;
            if (particle_types.count(name) > 0
                || dispersion_relations.count(name) > 0
                || confining_potentials.count(name) > 0
                || pairwise_potentials.count(name) > 0) {
                std::cerr << "ERROR: Command";
                for (const auto &word : words) { std::cerr << ' ' << word; }
                std::cerr << " redefines a name (" << name
                          << ") which has already been defined." << std::endl;
                std::exit(EXIT_FAILURE);
            }
            return name;
        }

    public: // =================================================================

        void add_particle(const ScriptCommand &command) {
            // =================================================== VALIDATE TYPE
            const std::string &type_str = get_string_parameter(command, "type");
            const auto type_iter = particle_types.find(type_str);
            if (type_iter == particle_types.end()) {
                std::cerr << "ERROR: Command " << command
                          << " references undeclared particle type \""
                          << type_str << "\"." << std::endl;
                std::exit(EXIT_FAILURE);
            }
            // Note: CLion does not yet understand C++17 structured bindings.
            // @formatter:off
            const auto [type_id, type_stat, type_spin] = type_iter->second;
            // @formatter:on
            if (type_stat.type != ExchangeStatistics::Type::FERMION) {
                std::cerr << "ERROR: Currently, only particles with "
                          << "fermionic exchange statistics are supported. "
                          << "Bosonic and distinguishable particle types "
                          << "will be supported in a future release of ZSVM."
                          << std::endl;
                std::exit(EXIT_FAILURE);
            }
            // =================================================== VALIDATE SPIN
            const auto spin_value = static_cast<int>(std::round(
                    2.0 * get_double_parameter(command, "spin")));
            if (std::abs(spin_value) > type_spin) {
                std::cerr << "ERROR: Cannot add particle of type \""
                          << type_str << "\" with z-component spin "
                          << to_signed_spin(spin_value)
                          << " exceeding the particle type's total spin "
                          << to_unsigned_spin(type_spin) << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (std::abs(spin_value) % 2 != type_spin % 2) {
                std::cerr << "ERROR: Cannot add particle of type \""
                          << type_str << "\" with z-component spin "
                          << to_signed_spin(spin_value)
                          << " of parity not matching its particle type's "
                          << "total spin " << to_unsigned_spin(type_spin)
                          << ". Both must be integers or half-integers."
                          << std::endl;
                std::exit(EXIT_FAILURE);
            }
            // ================================================ EXTRACT CARRIERS
            std::map<std::string, double> carriers;
            for (const auto &pair : command.get_named_parameters()) {
                if (pair.first != "type" && pair.first != "spin") {
                    if (pair.second.type != ScriptToken::Type::INTEGER
                        && pair.second.type != ScriptToken::Type::DECIMAL) {
                        std::cerr << "ERROR: Particle carrier \""
                                  << pair.first << "\" is expected to have "
                                  << "an integer or real value." << std::endl;
                        std::exit(EXIT_FAILURE);
                    }
                    carriers.insert({pair.first, pair.second.double_value});
                }
            }
            // =============================================== REGISTER PARTICLE
            particles.emplace_back(type_id, spin_value, carriers);
            // =================================================== PRINT SUMMARY
            std::cout << "INFO: Adding particle with type ID "
                      << type_id << ", "
                      << to_descriptive_string(type_stat.type)
                      << ", and z-component spin "
                      << to_signed_spin(spin_value) << '.' << std::endl;
            for (const auto &pair : carriers) {
                std::cout << "    " << pair.first
                          << " = " << pair.second << std::endl;
            }
        }

        void declare_particle_type(const ScriptCommand &command) {
            // Precondition: command has at least two words, and
            // the first two words are <DECLARE> <PARTICLE_TYPE>.
            const std::string &name_value = get_unique_identifier(command);
            const auto spin_value = static_cast<int>(std::round(
                    2.0 * get_double_parameter(command, "spin")));
            if (spin_value < 0) {
                std::cerr << "ERROR: Cannot declare particle type \""
                          << name_value << "\" with negative total spin."
                          << std::endl;
                std::exit(EXIT_FAILURE);
            }
            const ExchangeStatistics stat_value(get_enum_parameter(
                    command, "statistics", ExchangeStatistics::map()));
            const std::size_t id = particle_types.size();
            particle_types.insert({name_value,
                                   std::tie(id, stat_value, spin_value)});
            std::cout << "INFO: Registered particle type \""
                      << name_value << "\" with "
                      << to_descriptive_string(stat_value.type)
                      << " and spin " << to_unsigned_spin(spin_value)
                      << '.' << std::endl;
        }

        void declare_dispersion_relation(const ScriptCommand &command) {
            // Precondition: command has at least two words, and
            // the first two words are <DECLARE> <DISPERSION_RELATION>.
            const std::string &name = get_unique_identifier(command);
            const DispersionRelation<double>::Type type = get_enum_parameter(
                    command, "interaction", DispersionRelation<double>::map());
            const double strength = get_double_parameter(command, "strength");
            const double exponent = get_double_parameter(command, "exponent");
            const std::string &carrier =
                    get_string_parameter(command, "carrier");
            const DispersionRelation<double> dispersion_relation(
                    type, strength, exponent, carrier);
            dispersion_relations.insert({name, dispersion_relation});
            // TODO: Print dispersion relation type.
            std::cout << "INFO: Registered dispersion relation \"" << name
                      << "\" with strength " << strength
                      << ", exponent " << exponent
                      << ", and carrier \"" << carrier << "\"." << std::endl;
        }

        void declare_confining_potential(const ScriptCommand &command) {
            // Precondition: command has at least two words, and
            // the first two words are <DECLARE> <CONFINING_POTENTIAL>.
            const std::string &name = get_unique_identifier(command);
            const ConfiningPotential<double>::Type type = get_enum_parameter(
                    command, "interaction", ConfiningPotential<double>::map());
            const double strength = get_double_parameter(command, "strength");
            const double exponent = get_double_parameter(command, "exponent");
            const std::string &carrier =
                    get_string_parameter(command, "carrier");
            const ConfiningPotential<double> confining_potential(
                    type, strength, exponent, carrier);
            confining_potentials.insert({name, confining_potential});
            // TODO: Print confining potential type.
            std::cout << "INFO: Registered confining potential \"" << name
                      << "\" with strength " << strength
                      << ", exponent " << exponent
                      << ", and carrier \"" << carrier << "\"." << std::endl;
        }

        void declare_pairwise_potential(const ScriptCommand &command) {
            // Precondition: command has at least two words, and
            // the first two words are <DECLARE> <PAIRWISE_POTENTIAL>.
            const std::string &name = get_unique_identifier(command);
            const PairwisePotential<double>::Type type = get_enum_parameter(
                    command, "interaction", PairwisePotential<double>::map());
            const double strength = get_double_parameter(command, "strength");
            const double exponent = get_double_parameter(command, "exponent");
            const std::string &carrier =
                    get_string_parameter(command, "carrier");
            const PairwisePotential<double> pairwise_potential(
                    type, strength, exponent, carrier);
            pairwise_potentials.insert({name, pairwise_potential});
            // TODO: Print pairwise potential type.
            std::cout << "INFO: Registered pairwise potential \"" << name
                      << "\" with strength " << strength
                      << ", exponent " << exponent
                      << ", and carrier \"" << carrier << "\"." << std::endl;
        }

        bool is_jacobi_reducible() {
            if (dispersion_relations.size() != 1
                || !confining_potentials.empty()
                || pairwise_potentials.size() != 1) {
                return false;
            }
            const auto &dispersion = dispersion_relations.begin()->second;
            if (dispersion.type !=
                DispersionRelation<double>::Type::RADIAL_POWER_LAW
                || dispersion.strength != 0.5
                || dispersion.exponent != 2) {
                return false;
            }
            const auto &pairwise = pairwise_potentials.begin()->second;
            if (pairwise.type !=
                PairwisePotential<double>::Type::RADIAL_POWER_LAW
                || pairwise.exponent != -1
                || pairwise.strength != 1.0) {
                return false;
            }
            return true;
        }

        void expand_basis(const ScriptCommand &command) {
            // Precondition: command has at least two words, and
            // the first two words are <EXPAND> <BASIS>.
            if (!space_dimension.has_value()) {
                std::cerr << "ERROR: Space dimension must be set before "
                          << "issusing an \"expand basis\" command."
                          << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (dispersion_relations.empty()) {
                std::cerr << "WARNING: No dispersion relation terms have been "
                          << "declared. Physically nonsensical results may be "
                          << "produced." << std::endl;
            }
            if (confining_potentials.empty() && pairwise_potentials.empty()) {
                std::cerr << "WARNING: No potential energy terms have been "
                          << "declared. Physically nonsensical results may be "
                          << "produced." << std::endl;
            }
            const std::size_t target_size =
                    get_unsigned_parameter(command, "target_size");
            const std::size_t trials =
                    get_unsigned_parameter(command, "trials");
            const std::size_t max_iterations =
                    get_unsigned_parameter(command, "max_iterations");
            const std::size_t num_requested_threads =
                    get_unsigned_parameter(command, "num_threads");
            omp_set_num_threads(static_cast<int>(num_requested_threads));
            if (is_jacobi_reducible()) {
                std::cout << "Using Jacobi coordinates to reduce dimension.\n";

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

                std::chrono::high_resolution_clock::time_point start =
                        std::chrono::high_resolution_clock::now();
                std::vector<std::vector<double>> best_basis;
                std::vector<std::vector<std::vector<double>>>
                        thread_bases(num_threads);
                std::vector<double> thread_energies(num_threads);
                for (std::size_t i = 0; i < target_size; ++i) {
#pragma omp parallel
                    {
                        SphericalECGJacobiVariationalOptimizer optimizer(
                                *space_dimension, particles,
                                dispersion_relations.begin()->second.carrier,
                                pairwise_potentials.begin()->second.carrier);
                        optimizer.basis = best_basis;
                        optimizer.recompute_basis_matrices();
                        optimizer.recompute_solver_matrices();
                        optimizer.expand_amoeba(trials / num_threads,
                                                1.0, max_iterations);
                        thread_bases[omp_get_thread_num()] = optimizer.basis;
                        thread_energies[omp_get_thread_num()] =
                                optimizer.get_ground_state_energy();
                    }
                    std::size_t min_index = static_cast<std::size_t>(
                            std::min_element(thread_energies.begin(),
                                             thread_energies.end()) -
                            thread_energies.begin());
                    best_basis = thread_bases[min_index];
                    std::cout << i << '\t' << thread_energies[min_index]
                              << std::endl;
                }
                std::chrono::high_resolution_clock::time_point stop =
                        std::chrono::high_resolution_clock::now();
                std::cout << "Time elapsed (seconds): "
                          << std::chrono::duration_cast<
                                  std::chrono::milliseconds>(
                                  stop - start).count() / 1000.0 << std::endl;
            } else {
                SphericalECGVariationalOptimizer optimizer(
                        *space_dimension, particles, dispersion_relations,
                        confining_potentials, pairwise_potentials);
                for (std::size_t i = 0; i < target_size; ++i) {
                    optimizer.expand_amoeba(trials, 1.0, max_iterations);
                    std::cout << optimizer.get_ground_state_energy()
                              << std::endl;
                }
            }
        }

        void set_space_dimension(const ScriptCommand &command) {
            // Precondition: command has at least two words, and
            // the first two words are <SET> <SPACE_DIMENSION>.
            // TODO: Report error if unnecessary named parameters are given.
            space_dimension = get_integer(command);
            std::cout << "INFO: Setting ambient space dimension to "
                      << *space_dimension << "." << std::endl;
        }

        void run() {
            using T = ScriptToken::Type;
            while (true) {
                ScriptCommand command = parser.get_next_command();
                if (command.empty()) { break; }
                const auto &words = command.get_words();
                if (words[0].type == T::ADD && words.size() >= 2) {
                    switch (words[1].type) {
                        case T::PARTICLE:
                            add_particle(command);
                            break;
                        default:
                            std::cout << "WARNING: Unknown \"add\" command "
                                      << command << ". Ignoring this command."
                                      << std::endl;
                    }
                } else if (words[0].type == T::DECLARE && words.size() >= 2) {
                    switch (words[1].type) {
                        case T::PARTICLE_TYPE:
                            declare_particle_type(command);
                            break;
                        case T::DISPERSION_RELATION:
                            declare_dispersion_relation(command);
                            break;
                        case T::CONFINING_POTENTIAL:
                            declare_confining_potential(command);
                            break;
                        case T::PAIRWISE_POTENTIAL:
                            declare_pairwise_potential(command);
                            break;
                        default:
                            std::cout << "WARNING: Unknown \"declare\" command "
                                      << command << ". Ignoring this command."
                                      << std::endl;
                    }
                } else if (words[0].type == T::EXPAND && words.size() >= 2) {
                    switch (words[1].type) {
                        case T::BASIS:
                            expand_basis(command);
                            break;
                        default:
                            std::cout << "WARNING: Unknown \"expand\" command "
                                      << command << ". Ignoring this command."
                                      << std::endl;
                    }
                } else if (words[0].type == T::SET) {
                    if (words.size() >= 2
                        && words[1].type == T::SPACE_DIMENSION) {
                        set_space_dimension(command);
                    } else {
                        std::cout << "WARNING: Unknown \"set\" command "
                                  << command << ". Ignoring this command."
                                  << std::endl;
                    }
                } else {
                    std::cout << "WARNING: Unknown command " << command
                              << ". Ignoring this command." << std::endl;
                }
            }
        }

    }; // class ScriptInterpreter

} // namespace zsvm


int main(int argc, char **argv) {
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " SCRIPT_FILE" << std::endl;
        return 1;
    }
    std::string script_file_name(argv[1]);
    {
        std::ifstream script_file(script_file_name);
        if (!script_file.good()) {
            std::cout << "Error: could not open script file "
                      << script_file_name << "." << std::endl;
            return 2;
        }
    }
    std::cout << std::scientific;
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
    zsvm::ScriptInterpreter interpreter(script_file_name);
    interpreter.run();
    return 0;
}
