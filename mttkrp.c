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
#include "hacoo.h"
#include "matrix.h"
#include "common.hpp"
#include <omp.h>
#include <stdio.h>
#include <cblas.h>

/* Parallel MTTKRP */
matrix_t *mttkrp(struct hacoo_tensor *h, matrix_t **u, unsigned int n)
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

        //if (tid == 0) {
            //printf("Number of threads: %d\n", nthreads);
        //}

        partials[tid] = new_matrix(h->dims[n], fmax);
        matrix_t *local_res = partials[tid];

        int chunk = (h->nbuckets + nthreads - 1) / nthreads;
        int start = tid * chunk;
        int end = (start + chunk > h->nbuckets) ? h->nbuckets : start + chunk;

        unsigned int *idx = (unsigned int *) MALLOC(h->ndims * sizeof(unsigned int));
        double *rank_vec = (double *) MALLOC(fmax * sizeof(double));

        // Loop over assigned bucket vectors
        for (int i = start; i < end; i++) {
            bucket_vector *vec = &h->buckets[i];
            if (vec->size == 0)
                continue;

            for (size_t j = 0; j < vec->size; j++) {
                struct hacoo_bucket *cur = &vec->data[j];

                // Get full index array from compressed HaCOO format
                hacoo_extract_index(cur, h->ndims, idx);

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
                cblas_daxpy(fmax, 1.0, rank_vec, 1, local_res->data + idx[n] * fmax, 1);
            }
        }

        FREE(rank_vec); // Free thread-local buffer
        FREE(idx);
    }

    /* Time merge step */
    double t_start = omp_get_wtime();

    // Merge all thread-local results into the global result
    /* Hybrid */
    /*#pragma omp parallel for
    for (int i = 0; i < h->dims[n]; i++) {
        for (int f = 0; f < fmax; f++) {
            double sum = 0.0;
            for (int t = 0; t < num_threads; t++) {
                sum += partials[t]->vals[i][f];
            }
            res->vals[i][f] = sum;
        }
    }*/

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

    /*for collapse ver */
    /*pragma omp parallel for collapse(2)
    for (int i = 0; i < h->dims[n]; i++) {
        for (int f = 0; f < fmax; f++) {
            for (int t = 0; t < num_threads; t++) {
                res->vals[i][f] += partials[t]->vals[i][f];
            }
        }
    }*/

    double t_end = omp_get_wtime();
    printf("Merge time: %.6f seconds\n", t_end - t_start);

    for (int t = 0; t < num_threads; t++) {
        free_matrix(partials[t]);
    }

    FREE(partials);

    return res;
}

matrix_t *mttkrp_serial(struct hacoo_tensor *h, matrix_t **u, unsigned int n)
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

                hacoo_extract_index(cur, h->ndims, idx);

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

// function to test mttkrp
void mttkrp_test(struct hacoo_tensor *t)
{

  // Create factor matrices
  double a[] = {1, 3, 5, 2, 4, 6};

  double b[] = {1, 4, 7, 2, 5, 8, 3, 6, 9};

  double c[] = {1, 2, 3, 4, 5, 6};

  // make an array of 3 matrices
  int num_matrices = 3;

  matrix_t **u = (matrix_t **)MALLOC(sizeof(matrix_t *) * num_matrices);

  u[0] = array_to_matrix(a, 2, 3);
  u[1] = array_to_matrix(b, 3, 3);
  u[2] = array_to_matrix(c, 2, 3);

  matrix_t *m;

  for (int i = 0; i < num_matrices; i++)
  {
    m = mttkrp(t, u, i);
    printf("\nMode-%d MTTKRP: \n", i);
    print_matrix(m);
  }

  // free factor matrices
  for (int i = 0; i < num_matrices; i++)
  {
    free_matrix(u[i]);
  }
  FREE(u);

  // free m
  free_matrix(m);
}