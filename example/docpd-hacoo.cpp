/**
 * @file docpd-haccoo.cpp
 * @author Robert Lowe <rlowe8@utm.edu>
 * @brief Perform document cpd analysis using Haccoo.
 * @version 0.1
 * @date 2025-09-15
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include <iostream>
#include <fstream>
#include <map>
#include "tokenize.hpp"
#include "../hacoo.h"

int main(int argc, char **argv) 
{
    if(argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <ngram-size> <vocab_file> <doc_file>" << std::endl;
        return -1;
    }

    // Get ngram size
    int ngram_size = std::stoi(argv[1]);

    // Load the vocabulary
    std::ifstream vocab_ifs(argv[2]);
    if(!vocab_ifs) {
        std::cerr << "Error opening vocabulary file: " << argv[1] << std::endl;
        return -1;
    }
    auto vocab = load_vocab(vocab_ifs);
    vocab_ifs.close();

    //create the document tensor
    unsigned int dims[ngram_size];
    for(int i=0; i<ngram_size; ++i) {
        dims[i] = vocab.size();
    }
    struct hacoo_tensor *doc_tensor = hacoo_alloc(ngram_size, dims, 1024, 70);

    // Read the document file
    std::ifstream doc_ifs(argv[3]);
    if(!doc_ifs) {
        std::cerr << "Error opening document file: " << argv[3] << std::endl;
        hacoo_free(doc_tensor);
        return -1;
    }
    std::vector<std::string> tokens = tokenize(doc_ifs);
    doc_ifs.close();

    // insert n-grams into the tensor
    for(auto it=tokens.begin(); it != tokens.end() - ngram_size + 1; ++it) {
        unsigned int indices[ngram_size];
        for(int i=0; i<ngram_size; ++i) {
            indices[i] = vocab[it[i]];
        }
        hacoo_set(doc_tensor, indices, hacoo_get(doc_tensor, indices) + 1.0);
    }

    print_tensor(doc_tensor);
}