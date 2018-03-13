#include "ScriptParser.hpp"


zsvm::ScriptParser::ScriptParser(const std::string &script_file_name)
        : tokenizer(script_file_name) {}


zsvm::ScriptCommand zsvm::ScriptParser::get_next_command() {
    std::vector<ScriptToken> tokens;
    while (true) {
        ScriptToken token = tokenizer.get_next_token();
        if (token.get_type() == ScriptToken::Type::END_OF_FILE) {
            if (tokens.empty()) {
                return ScriptCommand();
            } else {
                throw std::invalid_argument(
                        "Script file ends with non-semicolon-"
                                "terminated command");
            }
        } else if (token.get_type() == ScriptToken::Type::SEMICOLON) {
            return ScriptCommand(tokens);
        } else {
            tokens.push_back(token);
        }
    }
}
