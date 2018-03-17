// C++ standard library headers
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>

// Project-specific headers
#include "EnumTypes.hpp"
#include "SphericalECGContext.hpp"
#include "RealVariationalSolver.hpp"
#include "AmoebaOptimizer.hpp"
#include "ScriptParser.hpp"

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

        const std::size_t num_particles;
        const std::size_t num_pairs;
        SphericalECGContext context;
        RealVariationalSolver solver;
        std::vector<std::vector<double>> basis;
        std::vector<std::vector<double>> basis_matrices;

    public: // ===================================================== CONSTRUCTOR

        explicit SphericalECGVariationalOptimizer(
                const std::vector<Particle> &particles, int space_dimension)
                : num_particles(particles.size()),
                  num_pairs(num_particles * (num_particles - 1) / 2),
                  context(SphericalECGContext::create(
                          particles, space_dimension)) {}

    public: // =================================================================

        double get_ground_state_energy() {
            return solver.get_eigenvalue(0);
        }

    public: // =================================================================

        double augmented_ground_state_energy(
                const double *__restrict__ new_basis_matrix) {
            if (solver.empty()) {
                double overlap_element, hamiltonian_element;
                context.compute_matrix_elements(
                        overlap_element, hamiltonian_element,
                        new_basis_matrix, new_basis_matrix);
                return hamiltonian_element / overlap_element;
            } else {
                const std::size_t basis_size = basis_matrices.size();
                std::vector<double> new_overlap_column(basis_size + 1);
                std::vector<double> new_hamiltonian_column(basis_size + 1);
                for (std::size_t i = 0; i < basis_size; ++i) {
                    context.compute_matrix_elements(
                            new_overlap_column[i], new_hamiltonian_column[i],
                            basis_matrices[i].data(), new_basis_matrix);
                }
                context.compute_matrix_elements(
                        new_overlap_column[basis_size],
                        new_hamiltonian_column[basis_size],
                        new_basis_matrix, new_basis_matrix);
                return solver.minimum_augmented_eigenvalue(
                        new_overlap_column.data(),
                        new_hamiltonian_column.data());
            }
        }

        void expand_random(std::size_t num_trials) {
            std::vector<double> new_basis_element(num_pairs);
            std::vector<double> new_basis_matrix(num_pairs);
            std::vector<double> best_basis_element(num_pairs);
            std::vector<double> best_basis_matrix(num_pairs);
            double best_energy = solver.empty()
                                 ? std::numeric_limits<double>::max()
                                 : solver.get_eigenvalue(0);
            for (std::size_t trial = 0; trial < num_trials; ++trial) {
                context.random_correlation_coefficients(
                        new_basis_element.data());
                context.gaussian_parameter_matrix(
                        new_basis_element.data(), new_basis_matrix.data());
                const double new_energy =
                        augmented_ground_state_energy(new_basis_matrix.data());
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

        void expand_amoeba(std::size_t num_trials,
                           double initial_step_size,
                           std::size_t max_steps) {
            std::vector<double> new_basis_element(num_pairs);
            std::vector<double> best_basis_element(num_pairs);
            double best_energy = solver.empty()
                                 ? std::numeric_limits<double>::max()
                                 : solver.get_eigenvalue(0);
            for (std::size_t trial = 0; trial < num_trials; ++trial) {
                context.random_correlation_coefficients(
                        new_basis_element.data());
                const double new_energy = refine_amoeba(
                        new_basis_element.data(), initial_step_size, max_steps);
                if (new_energy < best_energy) {
                    std::cout << "Improved energy to "
                              << new_energy << std::endl;
                    best_basis_element = new_basis_element;
                    best_energy = new_energy;
                }
            }
            basis.push_back(best_basis_element);
            std::vector<double> best_basis_matrix(num_pairs);
            context.gaussian_parameter_matrix(
                    best_basis_element.data(), best_basis_matrix.data());
            basis_matrices.push_back(best_basis_matrix);
            recompute_solver_matrices();
        }

    public: // =================================================================

        double refine_amoeba(double *__restrict__ basis_element,
                             double initial_step_size,
                             std::size_t max_steps) {
            std::vector<double> basis_matrix(num_pairs);
            dznl::AmoebaOptimizer amoeba(
                    basis_element, num_pairs, initial_step_size,
                    [&](const double *b) {
                        context.gaussian_parameter_matrix(
                                b, basis_matrix.data());
                        return augmented_ground_state_energy(
                                basis_matrix.data());
                    });
            for (std::size_t step = 0; step < max_steps; ++step) {
                amoeba.step();
            }
            return amoeba.current_minimum(basis_element);
        }

    public: // =================================================================

        void recompute_basis_matrices() {
            basis_matrices.clear();
            for (const auto &basis_element : basis) {
                std::vector<double> basis_matrix(num_pairs);
                context.gaussian_parameter_matrix(
                        basis_element.data(), basis_matrix.data());
                basis_matrices.push_back(basis_matrix);
            }
        }

        void recompute_solver_matrices() {
            solver.set_basis_size_destructive(basis_matrices.size());
            for (std::size_t i = 0; i < basis_matrices.size(); ++i) {
                for (std::size_t j = 0; j < basis_matrices.size(); ++j) {
                    context.compute_matrix_elements(
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


    struct GeneralParticle {

        const std::size_t type_id;
        const int spin;
        const std::map<std::string, double> carriers;

        GeneralParticle(std::size_t type_id, int spin,
                        const std::map<std::string, double> &carriers)
                : type_id(type_id), spin(spin), carriers(carriers) {}

    }; // class GeneralParticle


    class ScriptInterpreter {

        ScriptParser parser;
        std::optional<int> space_dimension;
        std::map<std::string, std::tuple<std::size_t, ExchangeStatistics, int>>
                particle_types; // name -> (ID, exchange statistics, spin)
        std::map<std::string, DispersionRelation> dispersion_relations;
        std::map<std::string, ConfiningPotential> confining_potentials;
        std::map<std::string, PairwisePotential> pairwise_potentials;
        std::vector<GeneralParticle> particles;

    public: // ===================================================== CONSTRUCTOR

        explicit ScriptInterpreter(const std::string &script_file_name)
                : parser(script_file_name),
                  space_dimension(),
                  particle_types(),
                  dispersion_relations(),
                  confining_potentials(),
                  pairwise_potentials(),
                  particles() {}

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
             // formatter:on
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
                    command, "statistics", ExchangeStatistics::MAP));
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
            const DispersionRelation::Type type = get_enum_parameter(
                    command, "interaction", DispersionRelation::MAP);
            const double strength = get_double_parameter(command, "strength");
            const double exponent = get_double_parameter(command, "exponent");
            const std::string &carrier =
                    get_string_parameter(command, "carrier");
            const DispersionRelation dispersion_relation(
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
            const ConfiningPotential::Type type = get_enum_parameter(
                    command, "interaction", ConfiningPotential::MAP);
            const double strength = get_double_parameter(command, "strength");
            const double exponent = get_double_parameter(command, "exponent");
            const std::string &carrier =
                    get_string_parameter(command, "carrier");
            const ConfiningPotential confining_potential(
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
            const PairwisePotential::Type type = get_enum_parameter(
                    command, "interaction", PairwisePotential::MAP);
            const double strength = get_double_parameter(command, "strength");
            const double exponent = get_double_parameter(command, "exponent");
            const std::string &carrier =
                    get_string_parameter(command, "carrier");
            const PairwisePotential pairwise_potential(
                    type, strength, exponent, carrier);
            pairwise_potentials.insert({name, pairwise_potential});
            // TODO: Print pairwise potential type.
            std::cout << "INFO: Registered pairwise potential \"" << name
                      << "\" with strength " << strength
                      << ", exponent " << exponent
                      << ", and carrier \"" << carrier << "\"." << std::endl;
        }

        void set_space_dimension(const ScriptCommand &command) {
            // Precondition: command has at least two words, and
            // the first two words are <SET> <SPACE_DIMENSION>.
            // TODO: Report error if unnecessary named parameters are given.
            // TODO: Make space_dimension a long long int.
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


int main() {
    std::cout << std::scientific;
    std::cout
            << std::setprecision(std::numeric_limits<double>::max_digits10);

    zsvm::ScriptInterpreter interpreter("../example_script.zscr");
    interpreter.run();

//    const zsvm::Particle electron_up = {0, 1.0, -1.0, zsvm::Spin::UP};
//    const zsvm::Particle electron_down = {0, 1.0, -1.0, zsvm::Spin::DOWN};
////    const zsvm::Particle positron_up = {1, 1.0, +1.0, zsvm::Spin::UP};
////    const zsvm::Particle positron_down = {1, 1.0, +1.0, zsvm::Spin::DOWN};
//    const zsvm::Particle beryllium_nucleus =
//            {2, 16538.028978017737, +4.0, zsvm::Spin::UP};
//    const std::vector<zsvm::Particle> particles = {
//            electron_up, electron_down, electron_up, electron_down,
//            beryllium_nucleus};
//    zsvm::SphericalECGVariationalOptimizer optimizer(particles, 3);
//
//    for (std::size_t basis_size = 0; basis_size < 10; ++basis_size) {
//        PRINT_EXECUTION_TIME(
//                optimizer.expand_amoeba(10, 0.5, 200);
//                optimizer.recompute_solver_matrices();
//                std::cout << optimizer.get_ground_state_energy() << std::endl;
//                std::ostringstream basis_output_file_name;
//                basis_output_file_name << "basis-";
//                basis_output_file_name << std::setw(8) << std::setfill('0');
//                basis_output_file_name << basis_size + 1 << ".tsv";
//                std::ofstream basis_output_file(basis_output_file_name.str());
//                basis_output_file << std::scientific;
//                basis_output_file << std::setprecision(
//                        std::numeric_limits<double>::max_digits10);
//                for (const auto &basis_element : optimizer.basis) {
//                    for (std::size_t i = 0; i < basis_element.size(); ++i) {
//                        if (i > 0) { basis_output_file << '\t'; }
//                        basis_output_file << basis_element[i];
//                    }
//                    basis_output_file << std::endl;
//                });
//    }
//
//#ifdef ZSVM_SPHERICAL_ECG_CONTEXT_TIMING_ENABLED
//    std::cout << "Matrix element calls:       "
//              << optimizer.context.get_matrix_element_calls() << std::endl;
//    std::cout << "Matrix element time:        "
//              << optimizer.context.get_matrix_element_time() << std::endl;
//#endif // ZSVM_SPHERICAL_ECG_CONTEXT_TIMING_ENABLED

    return 0;
}
