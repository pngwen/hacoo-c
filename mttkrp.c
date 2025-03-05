/*
Carry out MTTKRP between the tensor and an array of matrices,
unfolding the tensor along mode n.

Parameters:
  h - A pointer to a hacoo tensor with some nonnew_matrix.

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

matrix_t *mttkrp(struct hacoo_tensor *h, matrix_t **u, unsigned int n) {

    // Number of columns in factor matrices
    unsigned int fmax = u[0]->cols;

    // Create the global result array
    matrix_t *res = new_matrix(h->dims[n], fmax);

    // Gets the value from OMP_NUM_THREADS
    const int NUM_THREADS = omp_get_max_threads();
    printf("Max threads: %d\n", NUM_THREADS);
    
    omp_set_num_threads(NUM_THREADS);

    // Start OpenMP parallel region
    #pragma omp parallel
    {
        int num_threads = omp_get_num_threads();
        int thread_id = omp_get_thread_num();

        // Allocate thread-local variables
        unsigned int *idx = (unsigned int *)malloc(sizeof(unsigned int) * h->ndims);
        unsigned int *tind = (unsigned int *)malloc(sizeof(unsigned int) * h->nnz);
        double *t = (double *)malloc(sizeof(double) * h->nnz);
        matrix_t *local_res = new_matrix(h->dims[n], fmax); // Local result matrix

        if (!idx || !tind || !t || !local_res) {
            fprintf(stderr, "Error: Memory allocation failed.\n");
            free(idx);
            free(tind);
            free(t);
            free(local_res);
        }

        // Divide work among threads by columns
        for (int f = thread_id; f < fmax; f += num_threads) {
            int z = 0; // Counter for nonzeros

            for (int m = 0; m < h->nbuckets; m++) {
                if (!h->buckets[m]) {
                    continue; // Skip empty buckets
                }

                struct hacoo_bucket *cur;
                for (cur = h->buckets[m]; cur; cur = cur->next) {
                    hacoo_extract_index(cur, h->ndims, idx);

                    t[z] = cur->value;
                    tind[z] = idx[n];

                    for (int b = 0; b < h->ndims; b++) {
                        if (b == n) continue; // Skip the unfolded mode
                        t[z] *= u[b]->vals[idx[b]][f];
                    }

                    local_res->vals[tind[z]][f] += t[z];
                    z++;
                }
            }
        }

        // Merge local results into the global result
        #pragma omp critical
        {
            for (int i = 0; i < res->rows; i++) {
                for (int j = 0; j < res->cols; j++) {
                    res->vals[i][j] += local_res->vals[i][j];
                }
            }
        }

        // Free thread-local memory
        free(idx);
        free(tind);
        free(t);
        free(local_res);
    }

    return res;
}

matrix_t *mttkrp_serial(struct hacoo_tensor *h, matrix_t **u, unsigned int n) {

  // number of columns
  unsigned int fmax = u[0]->cols;

  // create the result array
  matrix_t *res = new_matrix(h->dims[n], fmax);

  // to hold current nnz's index
  unsigned int *idx = (unsigned int *)malloc(sizeof(unsigned int) * h->ndims);

  // tind holds index at specific mode f
  unsigned int *tind = (unsigned int *)malloc(sizeof(unsigned int) * h->nnz);

  // to hold nnz values
  double *t = (double *)malloc(sizeof(double) * h->nnz);

  //if (idx == NULL || tind == NULL || t == NULL) {
  if (tind == NULL || t == NULL) {
    fprintf(stderr, "Error: Memory allocation failed.\n");
    return NULL;
  }

  struct hacoo_bucket *cur, *next;

  // go through each column
  for (int f = 0; f < fmax; f++) {

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

        if (cur == NULL) {
          fprintf(stderr, "Error: cur is NULL.\n");
          return NULL;
        }

        if (z >= h->nnz) {
          fprintf(stderr,
                  "Error: z exceeds the number of non-zero elements.\n");
          return NULL;
        }

        t[z] = cur->value;
         if(t[z]==0) {printf("set t[%d] to 0.\n",z);}
        tind[z] = idx[n];

        int b = 0;
        
        while (b < fmax) {
          // skip the unfolded mode
          if (b == n) {
            b++;
            continue;
          }

          if (idx[b] >= u[b]->rows) {
            fprintf(stderr,
                    "Error: idx[%d] out of bounds for u[%d] with rows %d.\n", b,
                    b, u[b]->rows);
            return NULL;
          }

          // multiply the nnz by each factor
          t[z] *= u[b]->vals[idx[b]][f];
          b++; // advance to next factor matrix
        }
        z++; // advance to next nnz
      }
    }

    // accumulate m(:,f)
    for (int z = 0; z < h->nnz; z++) {
      if (tind[z] >= res->rows) {
        fprintf(stderr, "Error: tind[z] (%d) out of bounds for res with rows %d\n", tind[z], res->rows);
        return NULL;
      }

      res->vals[tind[z]][f] += t[z];
    }
  } // end for every column

  // free arrays
  free(idx);
  free(tind);
  free(t);

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

  u[0] = array_to_matrix(a, 2, 3);
  u[1] = array_to_matrix(b, 3, 3);
  u[2] = array_to_matrix(c, 2, 3);

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