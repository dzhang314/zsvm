#ifndef ZSVM_SCRIPT_TOKEN_HPP_INCLUDED
#define ZSVM_SCRIPT_TOKEN_HPP_INCLUDED

// C++ standard library headers
#include <iostream>
#include <map>
#include <string>

namespace zsvm {

    class ScriptToken {

    public: // =============================================== TOKEN ENUMERATION

        enum class Type {
            IDENTIFIER, INTEGER, DECIMAL,
            LEFT_PAREN, RIGHT_PAREN, EQUALS, COMMA, SEMICOLON,
            ADD, PARTICLE,
            DECLARE, PARTICLE_TYPE, DISPERSION_RELATION,
            CONFINING_POTENTIAL, PAIRWISE_POTENTIAL,
            EXPAND, AMOEBA, RANDOM,
            SET, SPACE_DIMENSION, BASIS_OUTPUT_FILE,
            ENERGY_OUTPUT_FILE, SUMMARY_OUTPUT_FILE,
            END_OF_FILE
        }; // enum class Type

        static const std::map<char, Type> PUNCTUATION_TOKENS;

        static const std::map<std::string, Type> KEYWORD_TOKENS;

    private: // =============================================== MEMBER VARIABLES

        Type type;
        std::string string_value;
        long long int integer_value;
        double double_value;

    public: // ===================================================== CONSTRUCTOR

        explicit ScriptToken(Type token_type);

        explicit ScriptToken(const std::string &identifier);

    public: // ======================================================= ACCESSORS

        Type get_type() const;

        std::string get_string_value() const;

        long long int get_integer_value() const;

        double get_double_value() const;

    }; // class ScriptToken

} // namespace zsvm

std::ostream &operator<<(std::ostream &os, const zsvm::ScriptToken &token);

#endif // ZSVM_SCRIPT_TOKEN_HPP_INCLUDED
