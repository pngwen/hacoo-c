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
#include <stdio.h>

matrix_t *mttkrp(struct hacoo_tensor *h, matrix_t **u, unsigned int n) {

  // number of columns
  unsigned int fmax = u[0]->cols;

  // create the result array
  matrix_t *res = zeroes(h->dims[n], fmax);

  // to hold current nnz's index
  unsigned int *idx = (unsigned int *)malloc(sizeof(unsigned int) * h->ndims);
  //printf("h->ndims: %d\n", h->ndims);
  
  // tind holds index at specific mode f
  unsigned int *tind = (unsigned int *)malloc(sizeof(unsigned int) * h->nnz);
  //printf("h->nnz: %d\n", h->nnz);
  // to hold nnz values
  double *t = (double *)malloc(sizeof(double) * h->nnz);

  if (idx == NULL || tind == NULL || t == NULL) {
    fprintf(stderr, "Error: Memory allocation failed.\n");
    return NULL;
  }

  struct hacoo_bucket *cur, *next;

  // go through each column
  for (int f = 0; f < fmax; f++) {
    //printf("fmax: %d\n", fmax);
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
        tind[z] = idx[n];

        int b = 0;
        /*
        printf("index: ");
        for (int k = 0; k < h->ndims; k++) {
          printf("%d ", idx[k]);
        }
        printf("\n");
        */
        while (b < fmax) {
          // skip the unfolded mode
          if (b == n) {
            b++;
            continue;
          }

          //printf("b: %d\n", b);
          if (idx[b] >= u[b]->rows) {
            fprintf(stderr,
                    "Error: idx[%d] out of bounds for u[%d] with rows %d.\n", b,
                    b, u[b]->rows);
            return NULL;
          }

          // multiply the nnz by each factor
          //printf("vals[idx[%d]][%d]: %f\n", z, f, u[b]->vals[idx[b]][f]);
          t[z] *= u[b]->vals[idx[b]][f];

          /*// print current t for debugging
          printf("t = ");
          for (int j = 0; j < h->ndims; j++) {
            printf("%f ", t[j]);
          }
          printf("\n");
          */

          b++; // advance to next factor matrix
        }
        z++; // advance to next nnz
        //printf("advancing to next nnz...\n");
      }
    }
    //printf("accumulate m(:,f)...\n");
    /*
    printf("tind:\n");
    //print tind
    for(int p = 0;p<h->nnz;p++){
      printf("%d ",tind[p]);
    }
    */
    
    // accumulate m(:,f)
    for (int z = 0; z < h->nnz; z++) {
      //printf("z: %d, tind[z]: %d, res->rows: %d\n", z, tind[z], res->rows);  // Debug: print indices and sizes
      if (tind[z] >= res->rows) {
        fprintf(stderr, "Error: tind[z] (%d) out of bounds for res with rows %d\n", tind[z], res->rows);
        return NULL;
      }
      
      res->vals[tind[z]][f] += t[z];
      // Debugging: print the updated result
      //printf("res[%d][%d] updated to: %f\n", tind[z], f, res->vals[tind[z]][f]);
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