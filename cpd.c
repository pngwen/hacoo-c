#include "cpd.h"
#include "hacoo.h"
#include "matrix.h"
#include "mttkrp.h"

//#define GRAMREG 1e-8
#define GRAMREG 0

void add_diagonal(matrix_t *matrix, double value) {
    for (unsigned int i = 0; i < matrix->rows && i < matrix->cols; i++) {
        matrix->vals[i][i] += value;
    }
}

// compute the gram product of a matrix
void gram_product(matrix_t *res, matrix_t **factor, unsigned int modes, unsigned int mode)
{
    matrix_t *g = new_matrix(factor[0]->cols, factor[0]->cols);
    matrix_t *res2 = new_matrix(res->rows, res->cols);

    // start with the identity matrix
    fill_identity_matrix(res2);

    for(int i=modes-1; i>=0; i--)
    {
        // skip the current mode
        if(i==mode) continue;
        fill_matrix(g, 0);
        mul_transpose_matrix(g, factor[i], factor[i]);
        fill_matrix(res, 0);
        mul_matrix(res, res2, g);
        for(int i=0; i<res->rows; i++)
        {
            for(int j=0; j<res->cols; j++)
            {
                res2->vals[i][j] = res->vals[i][j];
            }
        }
    }

    free_matrix(res2);
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
//        scale_matrix(factors[i], norm); 
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

//-- DEBUGGING
printf("Iter %u, mode %u: mttkrp_result norm = %f\n", iter, mode, matrix_frobenius_norm(mttkrp_result));
printf("Iter %u, mode %u: factor norm = %f\n", iter, mode, matrix_frobenius_norm(factors[mode]));
//-- END DEBUGGING

            free(mttkrp_result);
        }

        // Check for convergence (optional)
        // If converged, break the loop
    }

    free_matrix(gram);
    free_matrix(grami);

    return factors;
}
