#include "ScriptTokenizer.hpp"

// C++ standard library headers
#include <cmath> // for std::pow
#include <sstream> // for std::ostringstream
#include <stdexcept> // for std::invalid_argument
#include <vector> // for std::vector


zsvm::ScriptTokenizer::ScriptTokenizer(const std::string &script_file_name)
        : script_file(script_file_name) {}


bool zsvm::ScriptTokenizer::is_identifier_character(int c) {
    return std::isalnum(c) || c == '_';
}


bool zsvm::ScriptTokenizer::is_identifier_first_character(int c) {
    return std::isalpha(c) || c == '_';
}


bool zsvm::ScriptTokenizer::is_number_first_character(int c) {
    return std::isdigit(c) || c == '+' || c == '-';
}


zsvm::ScriptToken zsvm::ScriptTokenizer::get_next_token() {
    char next_char;
    if (!script_file.get(next_char)) {
        return ScriptToken(ScriptToken::Type::END_OF_FILE);
    }
    // At this point, we know we have successfully read a character.
    if (std::isspace(next_char)) {
        // Ignore consumed whitespace and try again.
        return get_next_token();
    }
    if (next_char == '#') {
        // Consume comment up to next line feed and try again.
        while (script_file.get(next_char)) {
            if (next_char == '\n') { break; }
        }
        return get_next_token();
    }
    // At this point, we have eliminated whitespace and comments,
    // and we expect to see the first character of a token.
    auto punct_iter = ScriptToken::PUNCTUATION_TOKENS.find(next_char);
    if (punct_iter != ScriptToken::PUNCTUATION_TOKENS.end()) {
        return ScriptToken(punct_iter->second);
    }
    if (is_identifier_first_character(next_char)) {
        std::vector<char> chars;
        chars.push_back(next_char);
        while (is_identifier_character(script_file.peek())) {
            script_file.get(next_char);
            chars.push_back(next_char);
        }
        std::string identifier(chars.begin(), chars.end());
        auto kw_iter = ScriptToken::KEYWORD_TOKENS.find(identifier);
        if (kw_iter == ScriptToken::KEYWORD_TOKENS.end()) {
            return ScriptToken(identifier);
        } else {
            return ScriptToken(kw_iter->second);
        }
    }
    if (is_number_first_character(next_char)) {
        const int sign = (next_char == '-') ? -1 : +1;
        long long int integer_value = std::isdigit(next_char)
                                      ? next_char - '0' : 0;
        while (std::isdigit(script_file.peek())) {
            script_file.get(next_char);
            integer_value = 10 * integer_value + (next_char - '0');
        }
        if (script_file.peek() == '.') {
            script_file.get(next_char); // Consume period.
            long long int fraction_value = 0, fraction_length = 0;
            while (std::isdigit(script_file.peek())) {
                script_file.get(next_char);
                fraction_value = 10 * fraction_value + (next_char - '0');
                ++fraction_length;
            }
            const double fraction = std::pow(10.0, -fraction_length)
                                    * fraction_value;
            const double decimal_value = integer_value + fraction;
            if (script_file.peek() == 'e' || script_file.peek() == 'E') {
                script_file.get(next_char); // Consume letter e.
                script_file.get(next_char); // Read next character.
                if (!is_number_first_character(next_char)) {
                    std::ostringstream error_message;
                    error_message << "Invalid character '" << next_char
                                  << "' in exponent of numeric constant";
                    throw std::invalid_argument(error_message.str());
                }
                const int exponent_sign = (next_char == '-') ? -1 : +1;
                if (next_char == '+' || next_char == '-') {
                    script_file.get(next_char);
                }
                if (!std::isdigit(next_char)) {
                    throw std::invalid_argument(
                            "Invalid exponent in numeric constant "
                                    "(expected to see a digit)");
                }
                long long int exponent_value = next_char - '0';
                while (std::isdigit(script_file.peek())) {
                    script_file.get(next_char);
                    exponent_value = 10 * exponent_value + (next_char - '0');
                }
                const double multiplier = std::pow(
                        10.0, exponent_sign * exponent_value);
                return ScriptToken(sign * decimal_value * multiplier);
            } else {
                return ScriptToken(sign * decimal_value);
            }
        } else if (script_file.peek() == 'e' || script_file.peek() == 'E') {
            script_file.get(next_char); // Consume letter e.
            script_file.get(next_char); // Read next character.
            if (!is_number_first_character(next_char)) {
                std::ostringstream error_message;
                error_message << "Invalid character '" << next_char
                              << "' in exponent of numeric constant";
                throw std::invalid_argument(error_message.str());
            }
            const int exponent_sign = (next_char == '-') ? -1 : +1;
            if (next_char == '+' || next_char == '-') {
                script_file.get(next_char);
            }
            if (!std::isdigit(next_char)) {
                throw std::invalid_argument(
                        "Invalid exponent in numeric constant "
                                "(expected to see a digit)");
            }
            long long int exponent_value = next_char - '0';
            while (std::isdigit(script_file.peek())) {
                script_file.get(next_char);
                exponent_value = 10 * exponent_value + (next_char - '0');
            }
            const double multiplier = std::pow(
                    10.0, exponent_sign * exponent_value);
            return ScriptToken(sign * integer_value * multiplier);
        } else {
            return ScriptToken(sign * integer_value);
        }
    }
    // At this point, we have exhausted all possible token types.
    std::ostringstream error_message;
    error_message << "Invalid character '" << next_char
                  << "' in script file";
    throw std::invalid_argument(error_message.str());
}
