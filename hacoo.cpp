/* File: hacoo.c
 * Purpose: Implementation of the hacoo sparse tensor library.
 */
#include "hacoo.h"
#include "alto.h"
#include "common.cpp"
#include "bitops.cpp"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include <immintrin.h>

#define LOAD 70
#define MIN_BUCKETS 128

/* Helper Function Prototypes */
static void hacoo_free_buckets(struct hacoo_tensor *t);

/* Calculate bucket index using packed ALTO encoding */
static size_t hacoo_bucket_index(struct hacoo_tensor *t, LIT packed);

/* Set hashing parameters */
static void hacoo_set_params(struct hacoo_tensor *t);

/* Search buckets if value already exists in the tensor*/
static struct hacoo_bucket *hacoo_search(bucket_vector *vec, LIT packed);
                                      
/* Return max number of bits used for encoding */
static size_t hacoo_max_bits(unsigned int n);

/* Allocation and deallocation functions */
struct hacoo_tensor *hacoo_alloc(unsigned int ndims, unsigned int *dims,
                                 size_t nbuckets, unsigned int load)
{
  struct hacoo_tensor *t = (struct hacoo_tensor *) malloc(sizeof(struct hacoo_tensor));

  if (t == NULL) {
    goto error;
  }

  /* initialize tensor fields */
  t->ndims = ndims;
  t->dims = (unsigned int *) MALLOC(ndims * sizeof(unsigned int));
  if (!t->dims) {
    goto error;
  }
  memcpy(t->dims, dims, sizeof(unsigned int) * ndims);

  t->nbuckets = nbuckets;
  t->load = load;
  t->nnz = 0;

  // Allocate array of bucket_vector structs (aligned for cache efficiency)
  t->buckets = (bucket_vector *)MALLOC(nbuckets * sizeof(bucket_vector));
  if (!t->buckets) {
      fprintf(stderr, "Error: Failed to allocate aligned memory for buckets.\n");
      goto error;
  }

  // Initialize each bucket_vector
  for (size_t i = 0; i < nbuckets; ++i) {
    t->buckets[i] = bucket_vector_create();
  }

  hacoo_set_params(t);

  //allocate mode masks
  t->mode_masks = (LIT*)calloc(t->ndims, sizeof(LIT));
  assert(t->mode_masks);
  
  //setup alto encoding
  alto_setup(t, LSB_FIRST, SHORT_FIRST);
  
  return t;

error:
  if (t) {
    hacoo_free(t);
  }
  return NULL;
}

void hacoo_free(struct hacoo_tensor *t) {
    if (!t) return;

    if (t->dims) {
        FREE(t->dims);
        t->dims = NULL;
    }
    if (t->buckets) {
        hacoo_free_buckets(t);
        t->buckets = NULL;
    }
    FREE(t);
}

/* Access functions */
void hacoo_set(struct hacoo_tensor *t, unsigned int *index, double value) {

  LIT alto_idx = alto_pack_index(index, t->mode_masks, t->ndims);
  size_t i = hacoo_bucket_index(t, alto_idx);

  bucket_vector *vec = &t->buckets[i];

  // Search for existing bucket with same packed index 
  struct hacoo_bucket *b = hacoo_search(vec, alto_idx);

  // If not found, insert new bucket
  if (!b) {
    struct hacoo_bucket new_bucket;
    new_bucket.alto_idx = alto_idx;
    new_bucket.value = value;

    bucket_vector_push_back(vec, new_bucket);
    t->nnz++; // Increment number of nonzeros
    return;
  }

  // If found, update value
  b->value = value;

  // Check if we need to rehash
  if (t->nbuckets > 0 &&
      ((double)t->nnz / (double)t->nbuckets) > ((double)t->load / 100.0)) {
    hacoo_rehash(&t);
    if (t == NULL) {
      fprintf(stderr, "Rehash failed, exiting.\n");
      return;
    }
  }
}

/* Set nonzero given alto_idx and value */
void hacoo_hset(struct hacoo_tensor *t, LIT alto_idx, double value) {

  size_t i = hacoo_bucket_index(t, alto_idx);

  bucket_vector *vec = &t->buckets[i];

  // Search for existing bucket with same packed index 
  struct hacoo_bucket *b = hacoo_search(vec, alto_idx);

  // If not found, insert new bucket
  if (!b) {
    struct hacoo_bucket new_bucket;
    new_bucket.alto_idx = alto_idx;
    new_bucket.value = value;

    bucket_vector_push_back(vec, new_bucket);
    t->nnz++; // Increment number of nonzeros
    return;
  }

  // If found, update value
  b->value = value;

  // Check if we need to rehash
  if (t->nbuckets > 0 &&
      ((double)t->nnz / (double)t->nbuckets) > ((double)t->load / 100.0)) {
    hacoo_rehash(&t);
    if (t == NULL) {
      fprintf(stderr, "Rehash failed, exiting.\n");
      return;
    }
  }
}

void hacoo_rehash(struct hacoo_tensor **t)
{
  // Step 1: Allocate new tensor with 2x buckets
  struct hacoo_tensor *dummy = hacoo_alloc((*t)->ndims, (*t)->dims, (*t)->nbuckets * 2, (*t)->load);
  if (dummy == NULL) {
    fprintf(stderr, "Failed to allocate dummy tensor during rehash.\n");
    return;
  }

  hacoo_set_params(dummy);

  unsigned int *index = (unsigned int *) malloc(sizeof(unsigned int) * (*t)->ndims);
  if (!index) {
    fprintf(stderr, "Error: Failed to allocate index array during rehash.\n");
    hacoo_free(dummy);
    return;
  }

  int nnz = 0;

  // Reinsert all elements from old tensor into new one
  for (size_t i = 0; i < (*t)->nbuckets; i++) {
    bucket_vector *vec = &(*t)->buckets[i];
    for (size_t j = 0; j < vec->size; ++j) {
      struct hacoo_bucket *b = &vec->data[j];
      alto_unpack(b->alto_idx, (*t)->mode_masks, (*t)->ndims, index);
      hacoo_set(dummy, index, b->value);
      nnz++;
    }
  }

  if ((*t)->nnz != nnz) {
    printf("Something went wrong. Only %d nnz copied when there were originally %d nnz.\n", nnz, (*t)->nnz);
  }

  // Copy important fields
  (*t)->sx = dummy->sx;
  (*t)->sy = dummy->sy;
  (*t)->sz = dummy->sz;

  // Free old bucket vectors
  hacoo_free_buckets(*t);

  // Swap buckets and metadata from dummy
  (*t)->buckets = dummy->buckets;
  (*t)->nbuckets = dummy->nbuckets;
  (*t)->nnz = dummy->nnz;

  dummy->buckets = NULL; // Prevent double free
  hacoo_free(dummy);
  free(index);
}

double hacoo_get(struct hacoo_tensor *t, unsigned int *index)
{
  LIT alto_idx = alto_pack_index(index, t->mode_masks, t->ndims);
  unsigned int i = hacoo_bucket_index(t, alto_idx);
  bucket_vector *vec = &t->buckets[i];

  // Search for existing bucket with same alto code
  struct hacoo_bucket *b = hacoo_search(vec, alto_idx);
  if (b)
  {
    return b->value;
  }

  return 0.0;
}

/* Helper function implementations. */

/* free buckets given a specific hacoo tensor*/
static void hacoo_free_buckets(struct hacoo_tensor *t)
{
  for (size_t i = 0; i < t->nbuckets; i++) {
    bucket_vector_free(&t->buckets[i]);
  }
  FREE(t->buckets);
  t->buckets = NULL;
}


static size_t hacoo_bucket_index(struct hacoo_tensor *t,
                                 LIT packed)
{
  unsigned long long hash = (unsigned long long) packed;

  hash = hash + (hash << t->sx);
  hash = hash ^ (hash >> t->sy);
  hash = hash + (hash << t->sz);
  return hash % t->nbuckets;
}

static void hacoo_set_params(struct hacoo_tensor *t)
{
  unsigned int bits;

  bits = ceil(log(t->nbuckets) / log(2));
  t->sx = ceil(bits / 8) - 1;
  t->sy = 4 * t->sx - 1;
  if (t->sy < 1)
  {
    t->sy = 1;
  }
  t->sz = ceil(bits / 2);
}


static struct hacoo_bucket *hacoo_search(bucket_vector *vec,
                                                LIT packed)
{
  for (size_t i = 0; i < vec->size; i++) {
    if (vec->data[i].alto_idx == packed) {
      return &vec->data[i];
    }
  }
  return NULL;
}

/* Return max number of bits used for encoding */
static size_t hacoo_max_bits(unsigned int n)
{
    size_t b1 = sizeof(uint64_t) * 8 / n;
    size_t b2 = sizeof(unsigned int) * 8;

    return b1 < b2 ? b1 : b2;
}

/*Allocate a new hacoo bucket.*/
struct hacoo_bucket *hacoo_new_bucket()
{
  struct hacoo_bucket *b;
  b = (struct hacoo_bucket *) MALLOC(sizeof(struct hacoo_bucket));
  if (!b) {
    fprintf(stderr, "Error: Failed to allocate memory for new bucket.\n");
    return NULL;
  }

  // Zero initialize
  memset(b, 0, sizeof(struct hacoo_bucket));

  return b;
}

/* Read the dimensions from stdin and build the tensor */
struct hacoo_tensor *hacoo_read_init()
{
  return hacoo_file_init(stdin);
}

/* Read an entry from stdin */
void hacoo_read_entry(struct hacoo_tensor *t)
{
  hacoo_file_entry(t, stdin);
}

/* Read a tensor from a tns file in COO format */
struct hacoo_tensor *hacoo_read_tensor_file(FILE *file)
{
  struct hacoo_tensor *t = hacoo_file_init(file);

  while(!feof(file)) {
    hacoo_file_entry(t, file);
  }

  return t;
}

/* Initialize a tensor from a file */
struct hacoo_tensor *hacoo_file_init(FILE *file) {

  // Buffer to read the input line
  char buffer[1024];
  fgets(buffer, sizeof(buffer), file);

  // Count the number of integers in the line
  unsigned int count = 0;
  for (char *p = buffer; *p; p++) {
    if (*p == ' ')
      count++;
  }
  count++;

  // Allocate memory for the array of dimensions
  unsigned int *dims = (unsigned int *) MALLOC(count * sizeof(unsigned int));
  if (!dims)
    return NULL;

  // Parse the input line and store integers in the array
  char *token = strtok(buffer, " ");
  for (unsigned int i = 0; i < count; i++) {
    dims[i] = strtoul(token, NULL, 10);
    token = strtok(NULL, " ");
  }

  // Allocate the tensor using the parsed dimensions
  struct hacoo_tensor *t = hacoo_alloc(count, dims, MIN_BUCKETS, LOAD);

  // Free the allocated memory for the dimensions array
  FREE(dims);

  return t;
}

/* Read an entry from a file */
void hacoo_file_entry(struct hacoo_tensor *t, FILE *file) {

  double value;
  unsigned int *index = (unsigned int *) malloc(t->ndims * sizeof(unsigned int));
  if (!index) {
    fprintf(stderr, "Error: Failed to allocate memory for index array.\n");
    return;
  }

  /* read the index- assumes indexes are 1-based */
  for (int i = 0; i < t->ndims; i++) {
    if (feof(file))
      return;
    fscanf(file, "%u", &index[i]);
    index[i] -= 1;
  }

  /* read the value */
  if (feof(file))
    return;
  fscanf(file, "%lf", &value);

  /* insert the value */
  hacoo_set(t, index, value);
}

/* Read a tensor file in HaCOO format */
struct hacoo_tensor *hacoo_read_htensor_file(FILE *file) {
  
  struct hacoo_tensor *t = hacoo_hfile_init(file);

  while(!feof(file)) {
    hacoo_hfile_entry(t, file);
  }

  return t;
}


/* Initialize a tensor from a file in HaCOO format */
struct hacoo_tensor *hacoo_hfile_init(FILE *file) {

    char buffer[1024];

    //read number of buckets
    if (fgets(buffer, sizeof(buffer), file) == NULL)
        return NULL;

    unsigned int numberOfBuckets = (unsigned int) strtoul(buffer, NULL, 10);

    //read dimensions line
    fgets(buffer, sizeof(buffer), file);

    // Count the number of integers in the line
    unsigned int count = 0;
    for (char *p = buffer; *p; p++) {
        if (*p == ' ')
        count++;
    }
    count++;

    // Allocate memory for the array of dimensions
    unsigned int *dims = (unsigned int *) MALLOC(count * sizeof(unsigned int));
    if (!dims)
        return NULL;

    // Parse the input line and store integers in the array
    char *token = strtok(buffer, " ");
    for (unsigned int i = 0; i < count; i++) {
        dims[i] = strtoul(token, NULL, 10);
        token = strtok(NULL, " ");
    } 

    struct hacoo_tensor *t = hacoo_alloc(count, dims, numberOfBuckets, LOAD);

    FREE(dims);

  return t;
}


/* Read an entry from a HaCOO format file */
void hacoo_hfile_entry(struct hacoo_tensor *t, FILE *file) {

  double value;
  LIT alto_idx = (LIT) malloc(sizeof(LIT));
  if (!alto_idx) {
    fprintf(stderr, "Error: Failed to allocate memory for ALTO index.\n");
    return;
  }

  /* read alto index */
  unsigned long tmp;
  fscanf(file, "%lu", &tmp);
  alto_idx = (LIT) tmp;

  /* read the value */
  if (feof(file))
    return;
  fscanf(file, "%lf", &value);

  /* insert the value */
  hacoo_hset(t, alto_idx, value);
}

/* Print out information about the tensor */
void hacoo_print_status(struct hacoo_tensor *t) {
  printf("nnz: %zu nbuckets: %zu\n", static_cast<long>(t->nnz), t->nbuckets);
}

void hacoo_print_tensor(struct hacoo_tensor *t)
{
    unsigned int index[t->ndims];       // for COO indices
    unsigned int coords[t->ndims];      // for unpacked ALTO indices

    for (size_t i = 0; i < t->nbuckets; i++) {
        bucket_vector *vec = &t->buckets[i];
        if (vec->size == 0) continue;

        printf("\nBucket %zu\n=============\n", i);

        for (size_t j = 0; j < vec->size; j++) {
            struct hacoo_bucket *b = &vec->data[j];

            // unpack ALTO index
            alto_unpack(b->alto_idx, t->mode_masks, t->ndims, coords);

            printf("0x%llu: ", static_cast <long long>(b->alto_idx));

            printf("COO indices: ");
            // print ALTO unpacked indices
            for (unsigned int k = 0; k < t->ndims; k++) {
                printf("%u ", coords[k]);
            }
            printf("\n");

            printf("%f\n", b->value);
        }
    }
}

/* Calculate the frobenius norm of the tensor */
double hacoo_frobenius_norm(struct hacoo_tensor *t)
{
    double norm = 0.0;
    for (size_t i = 0; i < t->nbuckets; i++) {
      bucket_vector *vec = &t->buckets[i];
      for (size_t j = 0; j < vec->size; j++) {
        struct hacoo_bucket *b = &vec->data[j];
        norm += b->value * b->value;
      }
    }
    return sqrt(norm);
}