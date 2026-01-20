/**
 * @file tokenize.hpp
 * @author Robert Lowe <rlowe8@utm.edu>
 * @brief Convert a string to a tokenized word.
 * @version 0.1
 * @date 2025-09-15
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include <cctype>
#include "tokenize.hpp"


// Tokenize a string by removing non-alphanumeric characters and converting to lowercase.
std::string tokenize(const std::string &input)
{
    std::string output;
    for (char c : input) {
        if (std::isalpha(c)) {
            output += std::tolower(c);
        }
    }
    return output;
}

// Tokenize an input stream into a list of tokens
std::vector<std::string> tokenize(std::istream &is)
{
    std::vector<std::string> tokens;
    std::string word;
    while (is >> word) {
        if(word.length() > 0) { 
            tokens.push_back(tokenize(word));
        }
    }
    return tokens;
}


// Load a vocabulary from a file stream. Format: word num
std::map<std::string, int> load_vocab(std::istream &is)
{
    std::map<std::string, int> vocab;
    std::string word;
    int index;
    while (is >> word >> index) {
        vocab[word] = index;
    }
    return vocab;
}