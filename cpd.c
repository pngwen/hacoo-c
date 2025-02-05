#include "cpd.h"
#include "hacoo.h"
#include "matrix.h"
#include "mttkrp.h"

// Function to solve the least squares problem AX = B
void least_squares(matrix_t *res, matrix_t *A, matrix_t *B)
{
    // Pre-allocated matrices for intermediate results
    matrix_t *A_transpose_A = new_matrix(A->cols, A->cols);
    matrix_t *A_transpose_B = new_matrix(A->cols, B->cols);

    // Compute A^T * A
    mul_transpose_matrix(A_transpose_A, A, A);

    // Compute the inverse of A^T * A
    invert_matrix(A_transpose_A, A_transpose_A);

    // Compute A^T * B
    mul_transpose_matrix(A_transpose_B, A, B);

    // Compute the solution X = (A^T * A)^(-1) * A^T * B
    mul_matrix(res, A_transpose_A, A_transpose_B);

    // Free intermediate matrices
    free_matrix(A_transpose_A);
    free_matrix(A_transpose_B);
}

// compute the canonical polyadic decomposition of a tensor
matrix_t **cpd(struct hacoo_tensor *t, unsigned int rank, unsigned int max_iter, double tol)
{
    // initialize matrices
    matrix_t **factors = malloc(t->ndims * sizeof(matrix_t *));

    for (unsigned int i = 0; i < t->ndims; i++)
    {
        factors[i] = new_random_matrix(t->dims[i], rank, 0, 1);
    }

    // solve the CPD via ALS
    for (unsigned int iter = 0; iter < max_iter; iter++)
    {
        for (unsigned int mode = 0; mode < t->ndims; mode++)
        {
            // Compute MTTKRP for the current mode
            matrix_t *mttkrp_result = mttkrp(t, factors, mode);

            // Solve the linear system to update the factor matrix
            least_squares(factors[mode], factors[mode], mttkrp_result);

            free(mttkrp_result);
        }

        // Check for convergence (optional)
        // If converged, break the loop
    }
}