#ifndef CPD_H
#define CPD_H
#include "hacoo.h"
#include "matrix.h"

/**
 * @brief Compute the canonical polyadic decomposition of a tensor.
 * 
 * @param t Pointer to the tensor to decompose
 * @param rank Number of factors to compute
 * @param max_iter Maximum number of iterations
 * @param tol Tolerance for convergence
 * @return matrix_t* List of rank matrices
 */
matrix_t **cpd(struct hacoo_tensor *t, unsigned int rank, unsigned int max_iter, double tol);

#endif