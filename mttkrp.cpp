/*
Carry out MTTKRP between the tensor and an array of matrices,
unfolding the tensor along mode n.

Parameters:
  h - A pointer to a hacoo tensor with some nonnew_matrix.

  u - A list of matrices, these correspond to the modes
    in the tensor, other than n. If i is the dimension in
    mode x, then u[x] must be an i x f matrix.
  n - The mode along which the tensor is unfolded for theac
    product.
Returns:
  A matrix with dimensions i_n x f
*/

#include "mttkrp.h"
#include "alto.h"
#include "hacoo.h"
#include "matrix.h"
#include "common.cpp"
#include <omp.h>
#include <cblas.h>
#include <stdio.h>
#include <cblas.h>

/* HaCOO Parallel MTTKRP */
matrix_t *hacoo_mttkrp(struct hacoo_tensor *h, matrix_t **u, unsigned int n)
{
    unsigned int fmax = u[0]->cols;

    // Allocate the final output matrix (global result)
    matrix_t *res = new_matrix(h->dims[n], fmax);

    int num_threads = omp_get_max_threads();

    matrix_t **partials = (matrix_t **) MALLOC(num_threads * sizeof(matrix_t *));

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int nthreads = omp_get_num_threads();
        int nnz_counter = 0;

        //if (tid == 0) {
            //printf("Number of threads: %d\n", nthreads);
        //}

        partials[tid] = new_matrix(h->dims[n], fmax);
        matrix_t *local_res = partials[tid];

        int chunk = (h->nbuckets + nthreads - 1) / nthreads;
        int start = tid * chunk;
        int end = (start + chunk > h->nbuckets) ? h->nbuckets : start + chunk;

        unsigned int *idx = (unsigned int *) malloc(h->ndims * sizeof(unsigned int));
        double *rank_vec = (double *) malloc(fmax * sizeof(double));

        // Loop over assigned bucket vectors
        for (int i = start; i < end; i++) {
            bucket_vector *vec = &h->buckets[i];
            if (vec->size == 0)
                continue;

            for (size_t j = 0; j < vec->size; j++) {
                nnz_counter++;
                struct hacoo_bucket *cur = &vec->data[j];

                // Get full index array from compressed HaCOO format

                // Get full index array from compressed ALTO format
                alto_unpack(cur->alto_idx, h->mode_masks, h->ndims, idx);

                // Initialize rank vector with cur->value
                for (int f = 0; f < fmax; f++) {
                    rank_vec[f] = cur->value;
                }

                // Multiply by the appropriate row from each factor matrix, skipping mode n
                for (int d = 0; d < h->ndims; d++) {
                    if (d == n) continue;
                    double *vec_d = u[d]->vals[idx[d]];
                    for (int f = 0; f < fmax; f++) {
                        rank_vec[f] *= vec_d[f];
                    }
                }

                // Accumulate into the local result row using daxpy
                for (int f = 0; f < fmax; f++) {
                    local_res->vals[idx[n]][f] += rank_vec[f];
                }
            }
        }

        free(rank_vec); // Free thread-local buffer
        free(idx);
    }

    // Merge all thread-local results into the global result
    /* Parallel over threads */
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int chunk = (h->dims[n] + num_threads - 1) / num_threads;
        int start = tid * chunk;
        int end = (start + chunk > h->dims[n]) ? h->dims[n] : start + chunk;

        for (int i = start; i < end; i++) {
            for (int f = 0; f < fmax; f++) {
                for (int t = 0; t < num_threads; t++) {
                    res->vals[i][f] += partials[t]->vals[i][f];
                }
            }
        }
    }

    for (int t = 0; t < num_threads; t++) {
        free_matrix(partials[t]);
    }

    FREE(partials);

    return res;
}

matrix_t *hacoo_mttkrp_serial(struct hacoo_tensor *h, matrix_t **u, unsigned int n)
{
    unsigned int fmax = u[0]->cols;
    matrix_t *res = new_matrix(h->dims[n], fmax);

    unsigned int *idx = (unsigned int *) MALLOC(sizeof(unsigned int) * h->ndims);
    unsigned int *tind = (unsigned int *) MALLOC(sizeof(unsigned int) * h->nnz);
    double *t = (double *) MALLOC(sizeof(double) * h->nnz);

    if (tind == NULL || t == NULL) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        return NULL;
    }

    for (int f = 0; f < fmax; f++) {
        int z = 0; // tracks the current nonzero

        for (int m = 0; m < h->nbuckets; m++) {
            bucket_vector *vec = &h->buckets[m];
            if (vec->size == 0) continue;

            for (size_t j = 0; j < vec->size; j++) {
                struct hacoo_bucket *cur = &vec->data[j];

                alto_unpack(cur->alto_idx, h->mode_masks, h->ndims, idx);


                if (cur == NULL) {
                    fprintf(stderr, "Error: cur is NULL.\n");
                    return NULL;
                }

                if (z >= h->nnz) {
                    fprintf(stderr, "Error: z exceeds nnz.\n");
                    return NULL;
                }

                t[z] = cur->value;
                tind[z] = idx[n];

                for (int d = 0; d < h->ndims; d++) {
                    if (d == n) continue;

                    if (idx[d] >= u[d]->rows) {
                        fprintf(stderr, "Error: idx[%d] out of bounds for u[%d] (rows = %d).\n",
                                d, d, u[d]->rows);
                        return NULL;
                    }

                    t[z] *= u[d]->vals[idx[d]][f];
                }

                z++;
            }
        }

        // Accumulate into output
        for (int z = 0; z < h->nnz; z++) {
            if (tind[z] >= res->rows) {
                fprintf(stderr, "Error: tind[%d] out of bounds for res (rows = %d).\n",
                        z, res->rows);
                return NULL;
            }

            res->vals[tind[z]][f] += t[z];
        }
    }

    FREE(idx);
    FREE(tind);
    FREE(t);

    return res;
}
