#include "cpd.h"
#include "hacoo.h"
#include "matrix.h"
#include "mttkrp.h"

#define GRAMREG 1e-8

void add_diagonal(matrix_t *matrix, double value) {
    for (unsigned int i = 0; i < matrix->rows && i < matrix->cols; i++) {
        matrix->vals[i][i] += value;
    }
}

// compute the gram product of a matrix
void gram_product(matrix_t *res, matrix_t **factor, unsigned int modes, unsigned int mode)
{
    matrix_t *g = new_matrix(factor[0]->cols, factor[0]->cols);

    // start with the identity matrix
    fill_identity_matrix(res);

    for(int i=modes-1; i>=0; i--)
    {
        // skip the current mode
        if(i==mode) continue;
        fill_matrix(g, 0);
        mul_transpose_matrix(g, factor[i], factor[i]);
        mul_matrix(res, res, g);
    }

    free_matrix(g);
}


// compute the canonical polyadic decomposition of a tensor
matrix_t **cpd(struct hacoo_tensor *t, unsigned int rank, unsigned int max_iter, double tol)
{
    // initialize matrices
    matrix_t **factors = malloc(t->ndims * sizeof(matrix_t *));
    matrix_t *gram = new_matrix(rank, rank);
    matrix_t *grami = new_matrix(rank, rank);
    double norm = frobenius_norm(t);

    for (unsigned int i = 0; i < t->ndims; i++)
    {
        factors[i] = new_random_matrix(t->dims[i], rank, 0, 1);
        scale_matrix(factors[i], norm); 
    }

    // solve the CPD via ALS
    for (unsigned int iter = 0; iter < max_iter; iter++)
    {
        for (unsigned int mode = 0; mode < t->ndims; mode++)
        {
            // Compute MTTKRP for the current mode
            matrix_t *mttkrp_result = mttkrp(t, factors, mode);

            // Compute the gram product and its inverse
            gram_product(gram, factors, t->ndims, mode);
            add_diagonal(gram, GRAMREG);
            invert_matrix(grami, gram);

            // Update the factor matrix
            fill_matrix(factors[mode], 0);
            mul_matrix(factors[mode], mttkrp_result, grami);

            free(mttkrp_result);
        }

        // Check for convergence (optional)
        // If converged, break the loop
    }

    free_matrix(gram);
    free_matrix(grami);

    return factors;
}
