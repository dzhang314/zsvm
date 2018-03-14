#include "ScriptCommand.hpp"

// C++ standard library headers
#include <stdexcept> // for std::invalid_argument


zsvm::ScriptToken zsvm::ScriptCommand::consume_first_token(
        std::vector<ScriptToken>::const_iterator &token_iterator,
        const std::vector<ScriptToken>::const_iterator &end_iterator) {
    if (token_iterator == end_iterator) {
        throw std::invalid_argument(
                "Command ended prematurely; cannot be empty");
    }
    const ScriptToken::Type type = token_iterator->type;
    if (type != ScriptToken::Type::ADD
        && type != ScriptToken::Type::DECLARE
        && type != ScriptToken::Type::EXPAND
        && type != ScriptToken::Type::SET) {
        throw std::invalid_argument(
                "Invalid command; first word must be add, declare, "
                        "expand, or set");
    }
    return *(token_iterator++);
}


std::string zsvm::ScriptCommand::consume_identifier(
        std::vector<ScriptToken>::const_iterator &token_iterator,
        const std::vector<ScriptToken>::const_iterator &end_iterator) {
    if (token_iterator == end_iterator) {
        throw std::invalid_argument(
                "Command ended prematurely; "
                        "expected to see parameter name");
    }
    if (token_iterator->type != ScriptToken::Type::IDENTIFIER) {
        throw std::invalid_argument(
                "Invalid command; "
                        "LHS of named parameter is not a valid name");
    }
    return (token_iterator++)->string_value;
}


void zsvm::ScriptCommand::consume_equals(
        std::vector<ScriptToken>::const_iterator &token_iterator,
        const std::vector<ScriptToken>::const_iterator &end_iterator) {
    if (token_iterator == end_iterator) {
        throw std::invalid_argument(
                "Command ended prematurely; "
                        "expected to see equals sign");
    }
    if (token_iterator->type != ScriptToken::Type::EQUALS) {
        throw std::invalid_argument(
                "Invalid command; LHS of named parameter "
                        "not followed by equals sign");
    }
    ++token_iterator;
}


zsvm::ScriptToken zsvm::ScriptCommand::consume_value(
        std::vector<ScriptToken>::const_iterator &token_iterator,
        const std::vector<ScriptToken>::const_iterator &end_iterator) {
    if (token_iterator == end_iterator) {
        throw std::invalid_argument(
                "Command ended prematurely; "
                        "expected to see named parameter value");
    }
    const ScriptToken::Type type = token_iterator->type;
    if (type != ScriptToken::Type::IDENTIFIER
        && type != ScriptToken::Type::INTEGER
        && type != ScriptToken::Type::DECIMAL) {
        throw std::invalid_argument(
                "Invalid command; "
                        "RHS of named parameter is not a valid value");
    }
    return *(token_iterator++);
}


bool zsvm::ScriptCommand::consume_right_paren_or_comma(
        std::vector<ScriptToken>::const_iterator &token_iterator,
        const std::vector<ScriptToken>::const_iterator &end_iterator) {
    if (token_iterator == end_iterator) {
        throw std::invalid_argument(
                "Command ended prematurely; "
                        "expected to see right parenthesis or comma");
    }
    const ScriptToken::Type type = token_iterator->type;
    if (type != ScriptToken::Type::RIGHT_PAREN
        && type != ScriptToken::Type::COMMA) {
        throw std::invalid_argument(
                "Invalid command; RHS of named parameter "
                        "not followed by right parenthesis or comma");
    }
    ++token_iterator;
    return (type == ScriptToken::Type::RIGHT_PAREN);
}


zsvm::ScriptCommand::ScriptCommand()
        : words(),
          named_parameters() {}


zsvm::ScriptCommand::ScriptCommand(const std::vector<ScriptToken> &tokens)
        : words(),
          named_parameters() {
    bool has_named_parameters = false;
    auto token_iterator = tokens.begin();
    words.push_back(consume_first_token(token_iterator, tokens.end()));
    for (; token_iterator != tokens.end(); ++token_iterator) {
        const ScriptToken::Type type = token_iterator->type;
        if (ScriptToken::PUNCTUATION_TOKENS.count(type)) {
            if (type == ScriptToken::Type::LEFT_PAREN) {
                ++token_iterator;
                has_named_parameters = true;
                break;
            } else {
                throw std::invalid_argument(
                        "Invalid command: "
                                "unexpected punctuation symbol");
            }
        } else {
            words.push_back(*token_iterator);
        }
    }
    if (has_named_parameters) {
        while (true) {
            std::string name = consume_identifier(
                    token_iterator, tokens.end());
            if (named_parameters.count(name)) {
                throw std::invalid_argument(
                        "Invalid command: duplicate named parameter");
            }
            consume_equals(token_iterator, tokens.end());
            ScriptToken value = consume_value(
                    token_iterator, tokens.end());
            named_parameters.insert({name, value});
            if (consume_right_paren_or_comma(
                    token_iterator, tokens.end())) {
                break;
            }
        }
        if (token_iterator != tokens.end()) {
            throw std::invalid_argument(
                    "Invalid command: continues after "
                            "closing right parenthesis");
        }
    }
}


bool zsvm::ScriptCommand::empty() const {
    // Note: a valid command must contain at least one word.
    return words.empty();
}


const std::vector<zsvm::ScriptToken> &zsvm::ScriptCommand::get_words() const {
    return words;
}


const std::map<std::string, zsvm::ScriptToken> &
zsvm::ScriptCommand::get_named_parameters() const {
    return named_parameters;
}


std::ostream &operator<<(std::ostream &os, const zsvm::ScriptCommand &command) {
    for (const auto &word : command.get_words()) {
        os << word << ' ';
    }
    os << '(';
    auto begin = command.get_named_parameters().begin();
    auto end = command.get_named_parameters().end();
    for (auto it = begin; it != end; ++it) {
        if (it != begin) { os << ", "; }
        os << it->first << " = " << it->second;
    }
    os << ')' << std::flush;
    return os;
}
