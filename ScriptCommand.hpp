#ifndef ZSVM_SCRIPT_COMMAND_HPP_INCLUDED
#define ZSVM_SCRIPT_COMMAND_HPP_INCLUDED

// C++ standard library headers
#include <map> // for std::map
#include <string> // for std::string
#include <vector> // for std::vector

// Project-specific headers
#include "ScriptToken.hpp"

namespace zsvm {

    class ScriptCommand {

    private: // =============================================== MEMBER VARIABLES

        std::vector<ScriptToken> words;
        std::map<std::string, ScriptToken> named_parameters;

    private: // ============================== STATIC CONSTRUCTOR HELPER METHODS

        static ScriptToken consume_first_token(
                std::vector<ScriptToken>::const_iterator &token_iterator,
                const std::vector<ScriptToken>::const_iterator &end_iterator);

        static std::string consume_identifier(
                std::vector<ScriptToken>::const_iterator &token_iterator,
                const std::vector<ScriptToken>::const_iterator &end_iterator);

        static void consume_equals(
                std::vector<ScriptToken>::const_iterator &token_iterator,
                const std::vector<ScriptToken>::const_iterator &end_iterator);

        static ScriptToken consume_value(
                std::vector<ScriptToken>::const_iterator &token_iterator,
                const std::vector<ScriptToken>::const_iterator &end_iterator);

        static bool consume_right_paren_or_comma(
                std::vector<ScriptToken>::const_iterator &token_iterator,
                const std::vector<ScriptToken>::const_iterator &end_iterator);

    public: // ===================================================== CONSTRUCTOR

        ScriptCommand();

        explicit ScriptCommand(const std::vector<ScriptToken> &tokens);

    public: // ======================================================= ACCESSORS

        bool empty() const;

        const std::vector<ScriptToken> &get_words() const;

        const std::map<std::string, ScriptToken> &get_named_parameters() const;

    }; // class ScriptCommand

} // namespace zsvm

std::ostream &operator<<(std::ostream &os, const zsvm::ScriptCommand &command);

#endif // ZSVM_SCRIPT_COMMAND_HPP_INCLUDED
