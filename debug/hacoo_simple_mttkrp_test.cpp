/* Basic test to test HaCOO MTTKRP. */

#include "hacoo.h"
#include "mttkrp.h"
#include <omp.h>
#include <cblas.h>
#include <stdio.h>
#include <stdlib.h>

void read_and_print(int argc, char *argv[]);

int main(int argc, char *argv[]) {
    // Force immediate output
    setvbuf(stdout, NULL, _IONBF, 0); // disable stdout buffering
    setvbuf(stderr, NULL, _IONBF, 0); // disable stderr buffering

    omp_set_num_threads(omp_get_max_threads()); 
    openblas_set_num_threads(omp_get_max_threads()); 

    FILE *file = fopen(argv[1], "r");

    if (argc < 2) {
        fprintf(stderr, "Usage: %s <tensor_file>\n", argv[0]);
        return -1;
    }
  
    if (!file) {
      perror("Error opening file");
      return -1;
    }
  
    // Read the tensor
    struct hacoo_tensor *t = hacoo_read_tensor_file(file);
    fclose(file);
  
    // Print tensor
    hacoo_print_tensor(t);

    hacoo_mttkrp_test(t);
    
    // Free tensor
    hacoo_free(t);

    return 0;
}


// function to test mttkrp
void hacoo_mttkrp_test(struct hacoo_tensor *t)
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

  print_matrix(u[0]);

  matrix_t *m;

  for (int i = 0; i < num_matrices; i++)
  {
    m = hacoo_mttkrp(t, u, i);
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