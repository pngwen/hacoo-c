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

/* Parallel MTTKRP */
matrix_t *hacoo_mttkrp(struct hacoo_tensor *h, matrix_t **u, unsigned int n)
{
    unsigned int fmax = u[0]->cols;

    // Allocate the final output matrix (global result)
    matrix_t *res = new_matrix(h->dims[n], fmax);

    int num_threads = omp_get_max_threads();
    matrix_t **partials = (matrix_t **) malloc(num_threads * sizeof(matrix_t *));

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int nthreads = omp_get_num_threads();

        partials[tid] = new_matrix(h->dims[n], fmax);
        matrix_t *local_res = partials[tid];

        int chunk = (h->nbuckets + nthreads - 1) / nthreads;
        int start = tid * chunk;
        int end = (start + chunk > h->nbuckets) ? h->nbuckets : start + chunk;

        unsigned int *idx = (unsigned int *) malloc(h->ndims * sizeof(unsigned int));

        for (int i = start; i < end; i++) {
            bucket_vector *vec = &h->buckets[i];
            if (vec->size == 0) continue;

            for (size_t j = 0; j < vec->size; j++) {
                struct hacoo_bucket *cur = &vec->data[j];

                // Unpack indices from ALTO encoding
                alto_unpack(cur->alto_idx, h->mode_masks, h->ndims, idx);

                double *out = local_res->vals[idx[n]];

                for (int f = 0; f < fmax; f++) {
                    double prod = cur->value;

                    for (int d = 0; d < h->ndims; d++) {
                        if (d == n) continue;
                        prod *= u[d]->vals[idx[d]][f];
                    }

                    out[f] += prod;
                }
            }
        }
        free(idx);
    }

    // Merge thread-local results into global matrix
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int chunk = (h->dims[n] + num_threads - 1) / num_threads;
        int start = tid * chunk;
        int end = (start + chunk > h->dims[n]) ? h->dims[n] : start + chunk;

        for (int i = start; i < end; i++) {
            for (int t = 0; t < num_threads; t++) {
                cblas_daxpy(fmax, 1.0, partials[t]->vals[i], 1, res->vals[i], 1);
            }
        }
    }

    for (int t = 0; t < num_threads; t++) {
        free_matrix(partials[t]);
    }
    free(partials);

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
