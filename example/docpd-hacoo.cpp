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
#include "hacoo.h"
#include "cpd.h"

void write_tsv_matrix(const std::string &filename, matrix_t *m);

int main(int argc, char **argv) 
{
    if(argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <ngram-size> <rank> <vocab_file> <doc_file>" << std::endl;
        return -1;
    }

    // Get ngram size
    int ngram_size = std::stoi(argv[1]);
    int rank = std::stoi(argv[2]);

    // Load the vocabulary
    std::ifstream vocab_ifs(argv[3]);
    if(!vocab_ifs) {
        std::cerr << "Error opening vocabulary file: " << argv[3] << std::endl;
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
    std::ifstream doc_ifs(argv[4]);
    if(!doc_ifs) {
        std::cerr << "Error opening document file: " << argv[4] << std::endl;
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

    // do the decomposition
    double tol = 1e-5;
    cpd_result_t *result= cpd(doc_tensor, rank, 1000, tol);

    // write the files
    for(unsigned int i=0; i < doc_tensor->ndims; ++i) {
        std::string filename = "factor_mode_" + std::to_string(i) + ".tsv";
        write_tsv_matrix(filename, result->factors[i]);
    }

    // write the lambda values
    {
        std::string filename = "lambdas.tsv";
        std::ofstream ofs(filename);
        if(ofs) {
            for(unsigned int i=0; i < result->rank; ++i) {
                ofs << result->lambda[i] << "\n";
            }
            ofs.close();
        } else {
            std::cerr << "Error opening output file: " << filename << std::endl;
        }
    }

    //clean up
    hacoo_free(doc_tensor);
    cpd_result_free(result);
}

void write_tsv_matrix(const std::string &filename, matrix_t *m)
{
    std::ofstream ofs(filename);
    if(!ofs) {
        std::cerr << "Error opening output file: " << filename << std::endl;
        return;
    }
    for(unsigned int i=0; i < m->rows; ++i) {
        for(unsigned int j=0; j < m->cols; ++j) {
            ofs << m->data[i*m->cols + j];
            if(j < m->cols - 1) {
                ofs << "\t";
            }
        }
        ofs << "\n";
    }
    ofs.close();
}