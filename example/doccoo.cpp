/**
 * @file docpd-haccoo.cpp
 * @author Robert Lowe <rlowe8@utm.edu>
 * @brief Convert a document to a COO document tensor.
 * @version 0.1
 * @date 2025-09-15
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <iterator>
#include <algorithm>
#include "tokenize.hpp"


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
        std::cerr << "Error opening vocabulary file: " << argv[3] << std::endl;
        return -1;
    }
    auto vocab = load_vocab(vocab_ifs);
    vocab_ifs.close();

    //create the document tensor
    std::vector<unsigned int> dims(ngram_size);
    for(int i=0; i<ngram_size; ++i) {
        dims[i] = vocab.size();
    }

    // Read the document file
    std::ifstream doc_ifs(argv[3]);
    if(!doc_ifs) {
        std::cerr << "Error opening document file: " << argv[3] << std::endl;
        return -1;
    }
    std::vector<std::string> tokens = tokenize(doc_ifs);
    doc_ifs.close();

    // insert n-grams into the tensor
    std::map<std::vector<unsigned int>, double> coo_data;
    for(auto it=tokens.begin(); it != tokens.end() - ngram_size + 1; ++it) {
        std::vector<unsigned int> indices(ngram_size);
        for(int i=0; i<ngram_size; ++i) {
            indices[i] = vocab[it[i]];
        }
        ++coo_data[indices];
    }

    // print the COO data
    std::ostream_iterator<unsigned int> out_it(std::cout, " ");
    std::copy(dims.begin(), dims.end()-1, out_it);
    std::cout << dims.back() << std::endl;
    for(const auto &entry : coo_data) {
        std::copy(entry.first.begin(), entry.first.end(), out_it);
        std::cout << entry.second << std::endl;
    }
}