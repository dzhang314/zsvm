#ifndef ZSVM_SCRIPT_TOKENIZER_HPP_INCLUDED
#define ZSVM_SCRIPT_TOKENIZER_HPP_INCLUDED

// C++ standard library headers
#include <fstream>

// Project-specific headers
#include "ScriptToken.hpp"

namespace zsvm {

    class ScriptTokenizer {

    private: // =============================================== MEMBER VARIABLES

        std::ifstream script_file;

    public: // ===================================================== CONSTRUCTOR

        explicit ScriptTokenizer(const std::string &script_file_name);

    private: // ========================================== STATIC HELPER METHODS

        static bool is_identifier_character(int c);

        static bool is_identifier_first_character(int c);

        static bool is_number_first_character(int c);

    public: // ======================================================= ACCESSORS

        ScriptToken get_next_token();

    }; // class ScriptTokenizer

} // namespace zsvm

#endif // ZSVM_SCRIPT_TOKENIZER_HPP_INCLUDED
