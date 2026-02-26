/* File: hacoo.h
 * Purpose: Header file for the hacoo sparse tensor library.
 */
#ifndef HACOO_H
#define HACOO_H
#include <stddef.h>
#include <stdint.h>
#include <stdint.h>
#include <stdio.h>
#include "vector.h"
#include "common.cpp"

struct hacoo_bucket {
  LIT alto_idx;  // packed ALTO encoding
  double value;
};

DEFINE_VECTOR_TYPE(struct hacoo_bucket, bucket_vector)

struct hacoo_tensor {
  size_t ndims;
  unsigned int *dims;
  bucket_vector *buckets; //vector of hacoo_buckets
  size_t nbuckets;
  unsigned int load;
  unsigned int nnz;
  unsigned int sx;
  unsigned int sy;
  unsigned int sz;
  LIT alto_mask;
  LIT *mode_masks;            // gather/scatter masks (nmode)
  #ifdef ALT_PEXT
    int *mode_pos;              // starting point for each mode mask (nmode)
  #endif
};

/* Allocation and deallocation functions */
struct hacoo_tensor *hacoo_alloc(unsigned int ndims, unsigned int *dims,
                                 size_t nbuckets, unsigned int load);
void hacoo_free(struct hacoo_tensor *t);

/* Rehash tensor that has exceeded load limit to new tensor */
void hacoo_rehash(struct hacoo_tensor **t);

/* Access functions */
void hacoo_set(struct hacoo_tensor *t, unsigned int *index, double value);
double hacoo_get(struct hacoo_tensor *t, unsigned int *index);

/* Allocate a new bucket */
struct hacoo_bucket *hacoo_new_bucket();

/* Read the dimensions from stdin and build the tensor */
struct hacoo_tensor *hacoo_read_init();

/* Read an entry from stdin */
void hacoo_read_entry(struct hacoo_tensor *t);

/* Read a tensor from a tns file in COO format */
struct hacoo_tensor *hacoo_read_tensor_file(FILE *file);

/* Initialize a tensor from a file in COO format */
struct hacoo_tensor *hacoo_file_init(FILE *file);

/* Read an entry from a file in COO format */
void hacoo_file_entry(struct hacoo_tensor *t, FILE *file);

/* Read tensor file in HaCOO format */
struct hacoo_tensor *hacoo_read_htensor_file(FILE *file);
struct hacoo_tensor *hacoo_hfile_init(FILE *file);
void hacoo_hfile_entry(struct hacoo_tensor *t, FILE *file);
void hacoo_hset(struct hacoo_tensor *t, LIT alto_idx, double value);

/* Print out information about the tensor */
void hacoo_print_status(struct hacoo_tensor *t);

/* Print the tensor hash table with COO listings */
void hacoo_print_tensor(struct hacoo_tensor *t);

/* Calculate the frobenius norm of the tensor */
double hacoo_frobenius_norm(struct hacoo_tensor *t);

#endif
