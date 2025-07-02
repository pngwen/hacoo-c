#include <math.h>
#include "cpd.h"
#include "hacoo.h"
#include "matrix.h"
#include "mttkrp.h"

#define GRAMREG 1e-6

// static helper prototypes
static void add_diagonal(matrix_t *matrix, double value);
static void gram_product(matrix_t *res, matrix_t **factor, unsigned int modes, unsigned int mode);
static cpd_result_t *cpd_alloc(struct hacoo_tensor *t, unsigned int rank);
static double normalize_column(matrix_t *m, unsigned int col_idx, unsigned int iter);
static void scale_factor_mode(cpd_result_t *result, unsigned int m, unsigned int iter);


static void add_diagonal(matrix_t *matrix, double value) {
    for (unsigned int i = 0; i < matrix->rows && i < matrix->cols; i++) {
        matrix->vals[i][i] += value;
    }
}

// compute the gram product of a matrix
static void gram_product(matrix_t *res, matrix_t **factor, unsigned int modes, unsigned int mode)
{
    matrix_t *g = new_matrix(factor[0]->cols, factor[0]->cols);

    // start with  matrix of ones
    for (unsigned int i = 0; i < res->rows; i++)
    {
        for (unsigned int j = 0; j < res->cols; j++)
        {
            res->vals[i][j] = 1.0;
        }
    }

    for(int n=0; n<modes; n++)
    {
        // skip the current mode
        if(n==mode) continue;
        mul_transpose_matrix(g, factor[n], factor[n]);

        // compute the hadamard product res .* g
        for (unsigned int j = 0; j < res->rows; j++)
        {
            for (unsigned int k = 0; k < res->cols; k++)
            {
                res->vals[j][k] *= g->vals[j][k];
            }
        }
    }

    free_matrix(g);
}


static cpd_result_t *cpd_alloc(struct hacoo_tensor *t, unsigned int rank)
{
    cpd_result_t *result = calloc(1, sizeof(cpd_result_t));
    if (!result) { goto bad; }

    result->rank = rank;
    result->ndims = t->ndims;

    // Allocate list of pointers for factor matrices
    result->factors = calloc(t->ndims, sizeof(matrix_t *));
    if(!result->factors) { goto bad; }

    // Allocate the lambda vector
    result->lambda = calloc(rank, sizeof(double));
    if (!result->lambda) { goto bad; }
    for(int i = 0; i < result->rank; i++)
    {
        result->lambda[i] = 1.0; // Initialize lambda to 1.0
    }

    // Initialize the random arrays
    for (unsigned int i = 0; i < t->ndims; i++)
    {
        result->factors[i] = new_random_matrix(t->dims[i], rank, 0, 1);
        if (!result->factors[i]) { goto bad; }
    }

    return result;

bad:
    cpd_result_free(result);
    return NULL;
}


/* Normalize a column of the matrix and return its l2 norm.*/
static double normalize_column(matrix_t *m, unsigned int col_idx, unsigned int iter)
{
    // compute the L2 norm
    double norm = 0.0;

    if(iter == 0) {
        // compute the L2 norm of the column
        for (unsigned int i = 0; i < m->rows; i++)
        {
            norm += m->vals[i][col_idx] * m->vals[i][col_idx];
        } 
        norm = sqrt(norm);
    } else {
        // compute the max-norm of the column
        norm = 1.0;
        for (unsigned int i = 0; i < m->rows; i++)
        {
            if (fabs(m->vals[i][col_idx]) > norm)
            {
                norm = fabs(m->vals[i][col_idx]);
            }
        }
    }

    // normalize the column
    for (unsigned int i = 0; i < m->rows; i++)
    {
        m->vals[i][col_idx] /= norm;
    }
    return norm;
}


static void scale_factor_mode(cpd_result_t *result, unsigned int m, unsigned int iter)
{
    for(unsigned int j=0; j<result->rank; j++){
        result->lambda[j] = normalize_column(result->factors[m], j, iter);
    }
}



// compute the canonical polyadic decomposition of a tensor
cpd_result_t *cpd(struct hacoo_tensor *t, unsigned int rank, unsigned int max_iter, double tol)
{
    // initialize matrices
    cpd_result_t *result = cpd_alloc(t, rank);
    matrix_t *gram = new_matrix(rank, rank);
    matrix_t *grami = new_matrix(rank, rank);
    double norm = frobenius_norm(t);


    // solve the CPD via ALS
    for (unsigned int iter = 0; iter < max_iter; iter++)
    {
        for (unsigned int mode = 0; mode < t->ndims; mode++)
        {
            // Compute MTTKRP for the current mode
            matrix_t *mttkrp_result = mttkrp(t, result->factors, mode);

            // Compute the gram product and its inverse
            gram_product(gram, result->factors, t->ndims, mode);
            //add_diagonal(gram, GRAMREG);
            invert_matrix(grami, gram);

            // Update the factor matrix
            mul_matrix(result->factors[mode], mttkrp_result, grami);
            scale_factor_mode(result, mode, iter);

//-- DEBUGGING
printf("Iter %u, mode %u: mttkrp_result norm = %f\n", iter, mode, matrix_frobenius_norm(mttkrp_result));
printf("Iter %u, mode %u: factor norm = %f\n", iter, mode, matrix_frobenius_norm(result->factors[mode]));
//-- END DEBUGGING

            free(mttkrp_result);
        }

        // Check for convergence (optional)
        // If converged, break the loop
    }

    free_matrix(gram);
    free_matrix(grami);

    return result;
}

// Free the memory allocated for the CPD result
void cpd_result_free(cpd_result_t *result)
{
    if(!result) return;

    if(result->factors) {
        for (unsigned int i = 0; i < result->ndims; i++) {
            free_matrix(result->factors[i]);
        }
        free(result->factors);
    }

    if(result->lambda) {
        free(result->lambda);
    }
    free(result);
}


void gram_test()
{
    matrix_t *a = new_matrix(2,2);
    matrix_t *g = new_matrix(2,2);

    /* populate a with:
       0.9027    0.4909
       0.9448    0.4893
     */
    a->vals[0][0] = 0.9027; a->vals[0][1] = 0.4909;
    a->vals[1][0] = 0.9448; a->vals[1][1] = 0.4893;

    printf("Iteration 1 Normalization:\n");
    copy_matrix_to(g,a);
    printf("lambda[0]=%f\n", normalize_column(g, 0, 0));
    printf("lambda[1]=%f\n", normalize_column(g, 1, 0));
    print_matrix(g);
    
    printf("Iteration 2 Normalization:\n");
    printf("lambda[0]=%f\n", normalize_column(a, 0, 1));
    printf("lambda[1]=%f\n", normalize_column(a, 1, 1));
    print_matrix(a);
}
