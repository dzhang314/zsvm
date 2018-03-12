#include "ScriptTokenizer.hpp"

// C++ standard library headers
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
        throw std::invalid_argument("NOT IMPLEMENTED");
    }
    // At this point, we have exhausted all possible token types.
    std::ostringstream error_message;
    error_message << "Invalid character '" << next_char
                  << "' in script file";
    throw std::invalid_argument(error_message.str());
}
