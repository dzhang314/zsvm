#ifndef ZSVM_SCRIPT_PARSER_HPP_INCLUDED
#define ZSVM_SCRIPT_PARSER_HPP_INCLUDED

// C++ standard library headers
#include <string>

// Project-specific headers
#include "ScriptCommand.hpp"
#include "ScriptTokenizer.hpp"

namespace zsvm {

    class ScriptParser {

        /* TODO: ScriptParser should eventually track line numbers and
         * positions for error reporting. */

    private: // =============================================== MEMBER VARIABLES

        ScriptTokenizer tokenizer;

    public: // ===================================================== CONSTRUCTOR

        explicit ScriptParser(const std::string &script_file_name);

    public: // ======================================================= ACCESSORS

        ScriptCommand get_next_command();

    }; // class ScriptParser

} // namespace zsvm

#endif // ZSVM_SCRIPT_PARSER_HPP_INCLUDED
