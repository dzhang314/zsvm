#include "ScriptTokenizer.hpp"

// C++ standard library headers
#include <cctype> // std::isalpha, std::isdigit, etc.
#include <cmath> // for std::pow
#include <sstream> // for std::ostringstream
#include <stdexcept> // for std::invalid_argument


zsvm::ScriptTokenizer::ScriptTokenizer(const std::string &script_file_name)
        : line_number(1),
          column_number(1),
          file_name(script_file_name),
          script_file(script_file_name) {}


bool zsvm::ScriptTokenizer::is_word_character(int c) {
    return std::isalnum(c) || c == '_';
}


bool zsvm::ScriptTokenizer::is_word_first_character(int c) {
    return std::isalpha(c) || c == '_';
}


bool zsvm::ScriptTokenizer::is_number_first_character(int c) {
    return std::isdigit(c) || c == '+' || c == '-';
}


bool zsvm::ScriptTokenizer::get(char &c) {
    script_file.get(c);
    if (c == '\n') {
        ++line_number;
        column_number = 1;
    } else {
        ++column_number;
    }
    return script_file.good();
}


zsvm::ScriptToken zsvm::ScriptTokenizer::read_word(char first) {
    std::vector<char> chars = {first};
    while (is_word_character(script_file.peek())) {
        char next;
        get(next);
        chars.push_back(next);
    }
    std::string identifier(chars.begin(), chars.end());
    auto kw_iter = ScriptToken::KEYWORD_TO_TOKEN_MAP.find(identifier);
    if (kw_iter == ScriptToken::KEYWORD_TO_TOKEN_MAP.end()) {
        return ScriptToken(line_number, column_number, identifier);
    } else {
        return ScriptToken(line_number, column_number, kw_iter->second);
    }
}


double zsvm::ScriptTokenizer::read_fraction(std::vector<char> &chars) {
    char next;
    get(next); // Consume period.
    chars.push_back(next);
    long long int fraction_value = 0, fraction_length = 0;
    while (std::isdigit(script_file.peek())) {
        get(next);
        chars.push_back(next);
        fraction_value = 10 * fraction_value + (next - '0');
        ++fraction_length;
    }
    return fraction_value / std::pow(10.0, fraction_length);
}


double zsvm::ScriptTokenizer::read_exponent(std::vector<char> &chars) {
    char next_char;
    get(next_char); // Consume letter e.
    chars.push_back(next_char);
    get(next_char); // Read next character.
    chars.push_back(next_char);
    if (!is_number_first_character(next_char)) {
        throw std::invalid_argument(
                "Invalid character in exponent of numeric constant");
    }
    const int exponent_sign = (next_char == '-') ? -1 : +1;
    if (next_char == '+' || next_char == '-') {
        get(next_char);
        chars.push_back(next_char);
    }
    if (!std::isdigit(next_char)) {
        throw std::invalid_argument(
                "Invalid exponent in numeric constant "
                        "(expected to see a digit)");
    }
    long long int exponent_value = next_char - '0';
    while (std::isdigit(script_file.peek())) {
        get(next_char);
        chars.push_back(next_char);
        exponent_value = 10 * exponent_value + (next_char - '0');
    }
    return std::pow(10.0, exponent_sign * exponent_value);
}


zsvm::ScriptToken zsvm::ScriptTokenizer::read_number(char first) {
    std::vector<char> chars = {first};
    char next_char;
    const int sign = (first == '-') ? -1 : +1;
    long long int integer_value = std::isdigit(first) ? first - '0' : 0;
    while (std::isdigit(script_file.peek())) {
        get(next_char);
        chars.push_back(next_char);
        integer_value = 10 * integer_value + (next_char - '0');
    }
    if (script_file.peek() == '.') {
        const double fraction = read_fraction(chars);
        const double decimal_value = integer_value + fraction;
        const double multiplier =
                (script_file.peek() == 'e' || script_file.peek() == 'E')
                ? read_exponent(chars) : 1.0;
        return ScriptToken(line_number, column_number,
                           sign * decimal_value * multiplier,
                           std::string(chars.begin(), chars.end()));
    } else if (script_file.peek() == 'e' || script_file.peek() == 'E') {
        const double multiplier = read_exponent(chars);
        return ScriptToken(line_number, column_number,
                           sign * integer_value * multiplier,
                           std::string(chars.begin(), chars.end()));
    } else {
        return ScriptToken(line_number, column_number,
                           sign * integer_value,
                           std::string(chars.begin(), chars.end()));
    }
}


zsvm::ScriptToken zsvm::ScriptTokenizer::get_next_token() {
    char next_char;
    if (!get(next_char)) {
        return ScriptToken(line_number, column_number,
                           ScriptToken::Type::END_OF_FILE);
    }
    // At this point, we know we have successfully read a character.
    // If we see a whitespace character, consume it and try again.
    if (std::isspace(next_char)) { return get_next_token(); }
    if (next_char == '#') {
        // Consume comment up to next line feed and try again.
        while (get(next_char)) { if (next_char == '\n') { break; }}
        return get_next_token();
    }
    // At this point, we have eliminated whitespace and comments,
    // and we expect to see the first character of a token.
    auto punct_iter = ScriptToken::PUNCTUATION_TO_TOKEN_MAP.find(next_char);
    if (punct_iter != ScriptToken::PUNCTUATION_TO_TOKEN_MAP.end()) {
        return ScriptToken(line_number, column_number, punct_iter->second);
    }
    if (is_word_first_character(next_char)) { return read_word(next_char); }
    if (is_number_first_character(next_char)) { return read_number(next_char); }
    // At this point, we have exhausted all possible token types.
    std::ostringstream error_message;
    error_message << "Invalid character '" << next_char << "' in script file";
    throw std::invalid_argument(error_message.str());
}
