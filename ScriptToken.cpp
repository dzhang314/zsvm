#include "ScriptToken.hpp"


const std::map<char, zsvm::ScriptToken::Type>
        zsvm::ScriptToken::PUNCTUATION_TO_TOKEN_MAP = {
        {'(', Type::LEFT_PAREN},
        {')', Type::RIGHT_PAREN},
        {'=', Type::EQUALS},
        {',', Type::COMMA},
        {';', Type::SEMICOLON},
};


const std::map<std::string, zsvm::ScriptToken::Type>
        zsvm::ScriptToken::KEYWORD_TO_TOKEN_MAP = {
        {"add",                 Type::ADD},
        {"particle",            Type::PARTICLE},
        {"declare",             Type::DECLARE},
        {"particle_type",       Type::PARTICLE_TYPE},
        {"dispersion_relation", Type::DISPERSION_RELATION},
        {"confining_potential", Type::CONFINING_POTENTIAL},
        {"pairwise_potential",  Type::PAIRWISE_POTENTIAL},
        {"expand",              Type::EXPAND},
        {"amoeba",              Type::AMOEBA},
        {"random",              Type::RANDOM},
        {"set",                 Type::SET},
        {"space_dimension",     Type::SPACE_DIMENSION},
        {"basis_output_file",   Type::BASIS_OUTPUT_FILE},
        {"energy_output_file",  Type::ENERGY_OUTPUT_FILE},
        {"summary_output_file", Type::SUMMARY_OUTPUT_FILE},
};


const std::set<zsvm::ScriptToken::Type>
        zsvm::ScriptToken::PUNCTUATION_TOKENS = {
        Type::LEFT_PAREN,
        Type::RIGHT_PAREN,
        Type::EQUALS,
        Type::COMMA,
        Type::SEMICOLON,
};


zsvm::ScriptToken::ScriptToken(
        std::size_t line_number, std::size_t column_number,
        Type token_type)
        : line_number(line_number),
          column_number(column_number),
          type(token_type),
          string_value(),
          integer_value(std::numeric_limits<long long int>::min()),
          double_value(std::numeric_limits<double>::quiet_NaN()) {}


zsvm::ScriptToken::ScriptToken(
        std::size_t line_number, std::size_t column_number,
        const std::string &identifier)
        : line_number(line_number),
          column_number(column_number),
          type(Type::IDENTIFIER),
          string_value(identifier),
          integer_value(std::numeric_limits<long long int>::min()),
          double_value(std::numeric_limits<double>::quiet_NaN()) {}


zsvm::ScriptToken::ScriptToken(
        std::size_t line_number, std::size_t column_number,
        long long int value, const std::string &repr)
        : line_number(line_number),
          column_number(column_number),
          type(Type::INTEGER),
          string_value(repr),
          integer_value(value),
          double_value(std::numeric_limits<double>::quiet_NaN()) {}


zsvm::ScriptToken::ScriptToken(
        std::size_t line_number, std::size_t column_number,
        double value, const std::string &repr)
        : line_number(line_number),
          column_number(column_number),
          type(Type::DECIMAL),
          string_value(repr),
          integer_value(std::numeric_limits<long long int>::min()),
          double_value(value) {}


std::ostream &operator<<(std::ostream &os, const zsvm::ScriptToken &token) {
    switch (token.type) {
        case zsvm::ScriptToken::Type::IDENTIFIER:
            os << "<IDENTIFIER, " << token.string_value << ">";
            break;
        case zsvm::ScriptToken::Type::INTEGER:
            os << "<INTEGER, " << token.integer_value << ">";
            break;
        case zsvm::ScriptToken::Type::DECIMAL:
            os << "<DECIMAL, " << token.double_value << ">";
            break;
        case zsvm::ScriptToken::Type::LEFT_PAREN:
            os << "<LEFT_PAREN>";
            break;
        case zsvm::ScriptToken::Type::RIGHT_PAREN:
            os << "<RIGHT_PAREN>";
            break;
        case zsvm::ScriptToken::Type::EQUALS:
            os << "<EQUALS>";
            break;
        case zsvm::ScriptToken::Type::COMMA:
            os << "<COMMA>";
            break;
        case zsvm::ScriptToken::Type::SEMICOLON:
            os << "<SEMICOLON>";
            break;
        case zsvm::ScriptToken::Type::ADD:
            os << "<ADD>";
            break;
        case zsvm::ScriptToken::Type::PARTICLE:
            os << "<PARTICLE>";
            break;
        case zsvm::ScriptToken::Type::DECLARE:
            os << "<DECLARE>";
            break;
        case zsvm::ScriptToken::Type::PARTICLE_TYPE:
            os << "<PARTICLE_TYPE>";
            break;
        case zsvm::ScriptToken::Type::DISPERSION_RELATION:
            os << "<DISPERSION_RELATION>";
            break;
        case zsvm::ScriptToken::Type::CONFINING_POTENTIAL:
            os << "<CONFINING_POTENTIAL>";
            break;
        case zsvm::ScriptToken::Type::PAIRWISE_POTENTIAL:
            os << "<PAIRWISE_POTENTIAL>";
            break;
        case zsvm::ScriptToken::Type::EXPAND:
            os << "<EXPAND>";
            break;
        case zsvm::ScriptToken::Type::AMOEBA:
            os << "<AMOEBA>";
            break;
        case zsvm::ScriptToken::Type::RANDOM:
            os << "<RANDOM>";
            break;
        case zsvm::ScriptToken::Type::SET:
            os << "<SET>";
            break;
        case zsvm::ScriptToken::Type::SPACE_DIMENSION:
            os << "<SPACE_DIMENSION>";
            break;
        case zsvm::ScriptToken::Type::BASIS_OUTPUT_FILE:
            os << "<BASIS_OUTPUT_FILE>";
            break;
        case zsvm::ScriptToken::Type::ENERGY_OUTPUT_FILE:
            os << "<ENERGY_OUTPUT_FILE>";
            break;
        case zsvm::ScriptToken::Type::SUMMARY_OUTPUT_FILE:
            os << "<SUMMARY_OUTPUT_FILE>";
            break;
        case zsvm::ScriptToken::Type::END_OF_FILE:
            os << "<END_OF_FILE>";
            break;
    }
    os << '(' << token.line_number << ", " << token.column_number << ')';
    return os;
}
