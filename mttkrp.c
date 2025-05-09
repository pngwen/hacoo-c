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
#include <omp.h>
#include <stdio.h>

matrix_t *mttkrp(struct hacoo_tensor *h, matrix_t **u, unsigned int n)
{
  unsigned int fmax = u[0]->cols;

  // Create the final output matrix (global result)
  matrix_t *res = new_matrix(h->dims[n], fmax);

  // Get the maximum number of threads available
  int num_threads = omp_get_max_threads();

  // Allocate an array of per-thread result matrices (no shared writes!)
  matrix_t **partials = malloc(num_threads * sizeof(matrix_t *));

// Start the parallel region
#pragma omp parallel
  {
    int tid = omp_get_thread_num();       // Thread ID
    int nthreads = omp_get_num_threads(); // Total threads

    // Partition the hash buckets among threads: [start, end)
    int chunk = (h->nbuckets + nthreads - 1) / nthreads;
    int start = tid * chunk;
    int end = (start + chunk > h->nbuckets) ? h->nbuckets : start + chunk;

    // Allocate a thread-local result matrix
    matrix_t *local_res = new_matrix(h->dims[n], fmax);
    partials[tid] = local_res;

    // Temporary buffer for tensor indices
    unsigned int *idx = malloc(h->ndims * sizeof(unsigned int));

    // Loop over the thread's assigned buckets
    for (int i = start; i < end; i++)
    {
      if (!h->buckets[i])
        continue;

      // Traverse the linked list of nonzeros in this bucket
      for (struct hacoo_bucket *cur = h->buckets[i]; cur; cur = cur->next)
      {
        hacoo_extract_index(cur, h->ndims, idx); // Extract indices of this nnz

        // Compute MTTKRP contribution for each column f
        for (int f = 0; f < fmax; f++)
        {
          double prod = cur->value; // Start with tensor value
          for (int d = 0; d < h->ndims; d++)
          {
            if (d == n)
              continue;                    // Skip the mode we're unfolding along
            prod *= u[d]->vals[idx[d]][f]; // Multiply by corresponding factor matrix entry
          }
          // Accumulate contribution into thread-local result matrix
          local_res->vals[idx[n]][f] += prod;
        }
      }
    }

    // Clean up thread-local index buffer
    free(idx);
  }

  // Merge per-thread results into the global result matrix
  #pragma omp parallel for collapse(2)
  for (int i = 0; i < h->dims[n]; i++) {
      for (int f = 0; f < fmax; f++) {
          for (int t = 0; t < num_threads; t++) {
              res->vals[i][f] += partials[t]->vals[i][f];
          }
      }
  }

  // Free the array of pointers to partial matrices
  free(partials);

  return res;
}

matrix_t *mttkrp_serial(struct hacoo_tensor *h, matrix_t **u, unsigned int n)
{

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

  // if (idx == NULL || tind == NULL || t == NULL) {
  if (tind == NULL || t == NULL)
  {
    fprintf(stderr, "Error: Memory allocation failed.\n");
    return NULL;
  }

  struct hacoo_bucket *cur, *next;

  // go through each column
  for (int f = 0; f < fmax; f++)
  {

    int z = 0; // counter for which nnz we're on

    for (int m = 0; m < h->nbuckets; m++)
    { // go through every bucket

      // if blank bucket, skip
      if (h->buckets[m] == NULL)
      {
        continue;
      }

      // go through each element in that bucket
      for (cur = h->buckets[m]; cur; cur = cur->next)
      {
        // decode element in the bucket
        hacoo_extract_index(cur, h->ndims, idx);

        if (cur == NULL)
        {
          fprintf(stderr, "Error: cur is NULL.\n");
          return NULL;
        }

        if (z >= h->nnz)
        {
          fprintf(stderr,
                  "Error: z exceeds the number of non-zero elements.\n");
          return NULL;
        }

        t[z] = cur->value;
        if (t[z] == 0)
        {
          printf("set t[%d] to 0.\n", z);
        }
        tind[z] = idx[n];

        int b = 0;

        while (b < fmax)
        {
          // skip the unfolded mode
          if (b == n)
          {
            b++;
            continue;
          }

          if (idx[b] >= u[b]->rows)
          {
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
    for (int z = 0; z < h->nnz; z++)
    {
      if (tind[z] >= res->rows)
      {
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
void mttkrp_test(struct hacoo_tensor *t)
{

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
  free(u);

  // free m
  free_matrix(m);
}

// for doubling the size of array
void resizeIntArray(unsigned int **arr, int originalSize)
{
  int newSize = originalSize * 2;
  int *temp = (int *)realloc(*arr, newSize * sizeof(int));

  if (temp == NULL)
  {
    printf("Memory reallocation failed.\n");
    exit(1);
  }

  *arr = temp;
}

void resizeDoubleArray(double **arr, int originalSize)
{
  int newSize = originalSize * 2;
  double *temp = (double *)realloc(*arr, newSize * sizeof(double));

  if (temp == NULL)
  {
    printf("Memory reallocation failed.\n");
    exit(1);
  }

  *arr = temp;
}