/**
 * @file vocab.cpp
 * @author Robert Lowe
 * @brief A simple vocabulary builder for text processing.
 * @version 0.1
 * @date 2025-09-15
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include <iostream>
#include <fstream>
#include <set>
#include <string>
#include "tokenize.hpp"

// Process the stream returning all unique lowercase words.
std::set<std::string> process_stream(std::istream &is);



int main(int argc, char **argv)
{
    std::set<std::string> vocab;

    // handle no arguments
    if(argc==1) {
        std::cerr << "Usage: " << argv[0] << " <file1> <file2> ..." << std::endl;
        return -1;
    }

    for(auto arg=argv+1; arg != argv+argc; ++arg) {
        std::ifstream ifs(*arg); if(!ifs) {
            std::cerr << "Error opening file: " << *arg << std::endl;
            continue;
        }
        auto file_vocab = process_stream(ifs);
        vocab.insert(file_vocab.begin(), file_vocab.end());
        ifs.close();
    }

    // print vocab
    int i=0;
    for(const auto &word : vocab) {
        std::cout << word << " " << ++i << std::endl;
    }
}


// Process the stream returning all unique lowercase words.
std::set<std::string> process_stream(std::istream &is)
{
    std::set<std::string> vocab;

    for (std::string &word : tokenize(is)) {
        if(word.length()>0) {
            vocab.insert(tokenize(word));
        }
    }

    return vocab;
}