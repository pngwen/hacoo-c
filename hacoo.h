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

struct hacoo_bucket {
  unsigned long long morton;
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
  //unsigned int base; //index base
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

/* extract the index from a bucket */
void hacoo_extract_index(struct hacoo_bucket *b, unsigned int n,
                         unsigned int *index);

/* Allocate a new bucket */
struct hacoo_bucket *hacoo_new_bucket();

/* Read the dimensions from stdin and build the tensor */
struct hacoo_tensor *read_init();

/* Read an entry from stdin */
void read_entry(struct hacoo_tensor *t);

/* Read a tensor from a tns file */
struct hacoo_tensor *read_tensor_file(FILE *file);

/* Delete this later */
struct hacoo_tensor *read_tensor_file_with_base(FILE *file, int zero_base);

/* Initialize a tensor from a file */
struct hacoo_tensor *file_init(FILE *file);

/* Read an entry from a file */
void file_entry(struct hacoo_tensor *t, FILE *file);

/* Delete this later*/
void file_entry_with_base(struct hacoo_tensor *t, FILE *file, int zero_base);

/* Print out information about the tensor */
void print_status(struct hacoo_tensor *t);

/* Print the tensor hash table with COO listings */
void print_tensor(struct hacoo_tensor *t);

/* Print the contents of a specific bucket in the tensor */
//void print_bucket(struct hacoo_tensor *t, int bucket_index);
//void print_bucket_from_ptr(struct hacoo_bucket *b, unsigned int ndims);
//void print_nth_nonzero(struct hacoo_tensor *t, int n);

/* Calculate the frobenius norm of the tensor */
double frobenius_norm(struct hacoo_tensor *t);

#endif
