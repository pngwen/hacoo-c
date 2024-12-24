/* File: hacoo.h
 * Purpose: Header file for the hacoo sparse tensor library.
 */
#ifndef HACOO_H
#define HACOO_H
#include <stddef.h>

struct hacoo_bucket {
  unsigned long long morton;
  double value;
  struct hacoo_bucket *next;
};

struct hacoo_tensor {
  unsigned int ndims;
  unsigned int *dims;
  struct hacoo_bucket **buckets;
  size_t nbuckets;
  unsigned int load;
  unsigned int nnz;
  unsigned int sx;
  unsigned int sy;
  unsigned int sz;
};

/* Allocation and deallocation functions */
struct hacoo_tensor *hacoo_alloc(unsigned int ndims, unsigned int *dims,
                                 size_t nbuckets, unsigned int load);
void hacoo_free(struct hacoo_tensor *t);

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

/* Print out information about the tensor */
void print_status(struct hacoo_tensor *t);

/* Print the tensor hash table with COO listings */
void print_tensor(struct hacoo_tensor *t);
#endif