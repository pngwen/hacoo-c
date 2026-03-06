/* Test to use ALTO encoding scheme */

#include "hacoo.h"
#include "alto.h"
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

    if (argc != 2) {
      printf("Usage: ./alto_encode_test <tensor_file>\n");
      return -1;
    }

    FILE *file = fopen(argv[1], "r");

    if (!file) {
      perror("Error opening file");
      return -1;
    }
  
    // Read the tensor
    struct hacoo_tensor *t = hacoo_read_tensor_file(file);
    fclose(file);
  
    // Print tensor status
    //print_status(t);

    // Print tensor
    hacoo_print_tensor(t);

    // Free tensor
    hacoo_free(t);

    return 0;
}