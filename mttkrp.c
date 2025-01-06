/*
Carry out MTTKRP between the tensor and an array of matrices,
unfolding the tensor along mode n.

Parameters:
  h - A pointer to a hacoo tensor with some nonzeroes.

  u - A list of matrices, these correspond to the modes
    in the tensor, other than n. If i is the dimension in
    mode x, then u[x] must be an i x f matrix.
  n - The mode along which the tensor is unfolded for the
    product.
Returns:
  A matrix with dimensions i_n x f
*/

#include "mttkrp.h"
#include "hacoo.h"
#include "matrix.h"
#include <omp.h>
#include <stdio.h>

#define NUM_THREADS 10

matrix_t *mttkrp(struct hacoo_tensor *h, matrix_t **u, unsigned int n) {

  omp_set_num_threads(NUM_THREADS);

  // number of columns
  unsigned int fmax = u[0]->cols;

  // create the global result array
  matrix_t *res = zeroes(h->dims[n], fmax);

// Parallelize the outer loop over factor matrix columns
#pragma omp parallel
  {

    int thread_id = omp_get_thread_num(); // Get the current thread ID

    if (thread_id == 0) {
      int num_threads = omp_get_num_threads();
      printf("Number of threads: %d\n", num_threads);
    }

// go through each column
#pragma omp for schedule(dynamic)
    for (int f = 0; f < fmax; f++) {
      // each thread gets its own cur, next pointers
      // every thread should also get its own idx, tind, and t
      // every thread should have a local result matrix to merge @ end
      // each thread will get a column to compute, so it doesn't have
      // to store the entire result matrix...

      // to hold current nnz's index
      unsigned int *idx =
          (unsigned int *)malloc(sizeof(unsigned int) * h->ndims);

      // tind holds index at specific mode f
      unsigned int *tind =
          (unsigned int *)malloc(sizeof(unsigned int) * h->nnz);

      // to hold nnz values
      double *t = (double *)malloc(sizeof(double) * h->nnz);

      // create local result array
      matrix_t *local_res = zeroes(h->dims[n], fmax);

      if (idx == NULL || tind == NULL || t == NULL || local_res == NULL) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        free(idx);
        free(tind);
        free(t);
        free(local_res);
        // return NULL;
      }

      struct hacoo_bucket *cur, *next;

      int z = 0; // counter for which nnz we're on

      for (int m = 0; m < h->nbuckets; m++) { // go through every bucket

        // if blank bucket, skip
        if (h->buckets[m] == NULL) {
          continue;
        }

        // go through each element in that bucket
        for (cur = h->buckets[m]; cur; cur = cur->next) {
          // decode element in the bucket
          hacoo_extract_index(cur, h->ndims, idx);

          t[z] = cur->value;
          tind[z] = idx[n];

          int b = 0;

          while (b < fmax) {
            // skip the unfolded mode
            if (b == n) {
              b++;
              continue;
            }

            // multiply the nnz by each factor
            t[z] *= u[b]->vals[idx[b]][f];

            b++; // advance to next factor matrix
          }
          // Accumulate into the local result matrix
          local_res->vals[tind[z]][f] += t[z];
          z++; // advance to next nnz
        }      // end for every element in that bucket
      }        // end for every bucket

      // accumulate m(:,f)
      // Reduce local result into global result
#pragma omp critical
      {
        for (int i = 0; i < res->rows; i++) {
          for (int j = 0; j < res->cols; j++) {
            res->vals[i][j] += local_res->vals[i][j];
          }
        }
      } // end critical
      // free arrays
      free(idx);
      free(tind);
      free(t);
      free(local_res);
    } // end for every column
  }   // end parallel

  return res;
}

// function to test mttkrp
void mttkrp_test(struct hacoo_tensor *t) {

  // Create factor matrices
  double a[] = {1, 3, 5, 2, 4, 6};

  double b[] = {1, 4, 7, 2, 5, 8, 3, 6, 9};

  double c[] = {1, 2, 3, 4, 5, 6};

  // make an array of 3 matrices
  int num_matrices = 3;

  matrix_t **u = (matrix_t **)malloc(sizeof(matrix_t *) * num_matrices);

  u[0] = new_matrix(a, 2, 3);
  u[1] = new_matrix(b, 3, 3);
  u[2] = new_matrix(c, 2, 3);

  matrix_t *m;

  for (int i = 0; i < num_matrices; i++) {
    m = mttkrp(t, u, i);
    printf("\nMode-%d MTTKRP: \n", i);
    print_matrix(m);
  }

  // free factor matrices
  for (int i = 0; i < num_matrices; i++) {
    free_matrix(u[i]);
  }
  free(u);

  // free m
  free_matrix(m);
}