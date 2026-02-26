/* Utility to convert COO file into HaCOO file
Format:
  number_of_buckets
  dim_1 dim_2 ... dim_n
  alto_encoding value
  alto_encoding value
  ...
  alto_encoding value
*/


#include "hacoo.h"
#include <omp.h>
#include <cblas.h>
#include <stdio.h>
#include <stdlib.h>

void read_and_print(int argc, char *argv[]);

int main(int argc, char *argv[]) {

  if (argc < 3) {
    fprintf(stderr, "Usage: %s <tensor_file> <output_file>\n", argv[0]);
    return -1;
  }

    FILE *infile = fopen(argv[1], "r");
  
    if (!infile) {
      perror("Error opening file");
      return -1;
    }
  
    // Read the tensor
    struct hacoo_tensor *t = hacoo_read_tensor_file(infile);
    fclose(infile);
  
    FILE *outfile = fopen(argv[2], "w");

    //write number of buckets
    fprintf(outfile,"%zu\n",t->nbuckets);

    //write dimensions
    for(int i=0;i<t->ndims;i++) {
      if(i == t->ndims-1) {
        fprintf(outfile, "%d", t->dims[i]);
      } else {
        fprintf(outfile, "%d ", t->dims[i]);
      }
    }
    fprintf(outfile,"\n");

    for (size_t i = 0; i < t->nbuckets; i++) {
      bucket_vector *vec = &t->buckets[i];
      //skip empty buckets
      if (vec->size == 0) continue;
      for (size_t j = 0; j < vec->size; j++) {
          struct hacoo_bucket *b = &vec->data[j];

          fprintf(outfile, "%llu %f\n",
                  (unsigned long long)b->alto_idx,
                  b->value);
      }
  }

  fclose(outfile);

  // Free tensor
  hacoo_free(t);

  return 0;
}