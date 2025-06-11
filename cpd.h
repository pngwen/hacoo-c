#ifndef CPD_H
#define CPD_H
#include "hacoo.h"
#include "matrix.h"

typedef struct cpd_result {
    unsigned int ndims;     // Number of modes of the tensor
    unsigned int rank;      // Rank of the decomposition
    matrix_t     **factors; // List of factor matrices
    double       *lambda;   // Scaling factors for each rank
} cpd_result_t;

/**
 * @brief Compute the canonical polyadic decomposition of a tensor.
 * 
 * @param t Pointer to the tensor to decompose
 * @param rank Number of factors to compute
 * @param max_iter Maximum number of iterations
 * @param tol Tolerance for convergence
 * @return matrix_t* List of rank matrices
 */
cpd_result_t *cpd(struct hacoo_tensor *t, unsigned int rank, unsigned int max_iter, double tol);

/**
 * @brief Free the memory allocated for the CPD result.
 * @param result Pointer to the cpd_result_t structure to free
 */
void cpd_result_free(cpd_result_t *result);
#endif
