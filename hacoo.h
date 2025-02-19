/* File: hacoo.h
 * Purpose: Header file for the hacoo sparse tensor library.
 */
#ifndef HACOO_H
#define HACOO_H
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

struct hacoo_bucket {
  unsigned long long morton;
  double value;
  struct hacoo_bucket *next;
};

struct hacoo_tensor {
  size_t ndims;
  unsigned int *dims;
  struct hacoo_bucket **buckets;
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

/* Initialize a tensor from a file */
struct hacoo_tensor *file_init(FILE *file);

/* Read an entry from a file */
void file_entry(struct hacoo_tensor *t, FILE *file);

/* Print out information about the tensor */
void print_status(struct hacoo_tensor *t);

/* Print the tensor hash table with COO listings */
void print_tensor(struct hacoo_tensor *t);

/* Print the contents of a specific bucket in the tensor */
void print_bucket(struct hacoo_tensor *t, int bucket_index);
void print_bucket_from_ptr(struct hacoo_bucket *b, unsigned int ndims);

void print_nth_nonzero(struct hacoo_tensor *t, int n);

#endif