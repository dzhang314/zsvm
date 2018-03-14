#ifndef ZSVM_SCRIPT_TOKENIZER_HPP_INCLUDED
#define ZSVM_SCRIPT_TOKENIZER_HPP_INCLUDED

// C++ standard library headers
#include <fstream> // for std::ifstream
#include <vector> // for std::vector

// Project-specific headers
#include "ScriptToken.hpp"

namespace zsvm {

    class ScriptTokenizer {

    private: // =============================================== MEMBER VARIABLES

        std::size_t line_number;
        std::size_t column_number;
        const std::string file_name;
        std::ifstream script_file;

    public: // ===================================================== CONSTRUCTOR

        explicit ScriptTokenizer(const std::string &script_file_name);

    private: // ================================================= HELPER METHODS

        static bool is_word_character(int c);

        static bool is_word_first_character(int c);

        static bool is_number_first_character(int c);

        bool get(char &c);

        ScriptToken read_word(char first);

        double read_fraction(std::vector<char> &chars);

        double read_exponent(std::vector<char> &chars);

        ScriptToken read_number(char first);

    public: // ======================================================= ACCESSORS

        ScriptToken get_next_token();

    }; // class ScriptTokenizer

} // namespace zsvm

#endif // ZSVM_SCRIPT_TOKENIZER_HPP_INCLUDED
