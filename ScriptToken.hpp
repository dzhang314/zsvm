#ifndef ZSVM_SCRIPT_TOKEN_HPP_INCLUDED
#define ZSVM_SCRIPT_TOKEN_HPP_INCLUDED

// C++ standard library headers
#include <iostream> // for std::ostream
#include <map> // for std::map
#include <set> // for std::set
#include <string> // for std::string

namespace zsvm {

    struct ScriptToken {

    public: // =============================================== TOKEN ENUMERATION

        enum class Type {
            IDENTIFIER, INTEGER, DECIMAL,
            LEFT_PAREN, RIGHT_PAREN, EQUALS, COMMA, SEMICOLON,
            ADD, PARTICLE,
            DECLARE, PARTICLE_TYPE, DISPERSION_RELATION,
            CONFINING_POTENTIAL, PAIRWISE_POTENTIAL,
            EXPAND, REFINE, CLEAR, BASIS,
            SET, SPACE_DIMENSION,
            BASIS_INPUT_FILE, BASIS_OUTPUT_FILE,
            ENERGY_OUTPUT_FILE, SUMMARY_OUTPUT_FILE,
            END_OF_FILE
        }; // enum class Type

        static const std::map<char, Type> PUNCTUATION_TO_TOKEN_MAP;

        static const std::map<std::string, Type> KEYWORD_TO_TOKEN_MAP;

        static const std::set<Type> PUNCTUATION_TOKENS;

    public: // ================================================ MEMBER VARIABLES

        const std::size_t line_number;
        const std::size_t column_number;
        const Type type;
        const std::string string_value;
        const long long int integer_value;
        const double double_value;

    public: // ===================================================== CONSTRUCTOR

        explicit ScriptToken(
                std::size_t line_number, std::size_t column_number,
                Type token_type);

        explicit ScriptToken(
                std::size_t line_number, std::size_t column_number,
                const std::string &identifier);

        explicit ScriptToken(
                std::size_t line_number, std::size_t column_number,
                long long int value, const std::string &repr);

        explicit ScriptToken(
                std::size_t line_number, std::size_t column_number,
                double value, const std::string &repr);

    }; // class ScriptToken

} // namespace zsvm

std::ostream &operator<<(std::ostream &os, const zsvm::ScriptToken &token);

#endif // ZSVM_SCRIPT_TOKEN_HPP_INCLUDED
