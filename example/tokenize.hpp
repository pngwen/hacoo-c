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
#pragma once
#include <string>
#include <map>
#include <iostream>
#include <vector>

/**
 * @brief Tokenize a string by removing non-alphanumeric characters and converting to lowercase.
 * 
 * @param input 
 * @return std::string 
 */
std::string tokenize(const std::string &input);

/**
 * @brief Tokenize an input stream into a list of tokens
 * 
 * @param is 
 * @return * std::vector<std::string> 
 */
std::vector<std::string> tokenize(std::istream &is);

/**
 * @brief Load a vocabulary from a file stream. Format: word num
 * 
 * @param is 
 * @return std::map<std::string, int> 
 */
std::map<std::string, int> load_vocab(std::istream &is);