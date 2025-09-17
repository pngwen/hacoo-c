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
static void free_buckets(bucket_vector *buckets, size_t nbuckets);
static size_t hacoo_bucket_index(struct hacoo_tensor *t,
                                 LIT packed);
static void hacoo_compute_params(struct hacoo_tensor *t);
static struct hacoo_bucket *hacoo_bucket_search(bucket_vector *vec,
                                                LIT packed);
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

  hacoo_compute_params(t);

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

void hacoo_free(struct hacoo_tensor *t)
{
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
void hacoo_set(struct hacoo_tensor *t, unsigned int *index, double value)
{

  LIT alto_idx = alto_pack_index(index, t->mode_masks, t->ndims);
  size_t i = hacoo_bucket_index(t, alto_idx);

  bucket_vector *vec = &t->buckets[i];

  // Search for existing bucket with same packed index 
  struct hacoo_bucket *b = hacoo_bucket_search(vec, alto_idx);

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

  hacoo_compute_params(dummy);

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
  struct hacoo_bucket *b = hacoo_bucket_search(vec, alto_idx);
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

/* free only buckets */
static void free_buckets(bucket_vector *buckets, size_t nbuckets)
{
  for (size_t i = 0; i < nbuckets; ++i) {
    bucket_vector_free(&buckets[i]);
  }
  FREE(buckets);
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

static void hacoo_compute_params(struct hacoo_tensor *t)
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


static struct hacoo_bucket *hacoo_bucket_search(bucket_vector *vec,
                                                LIT packed)
{
  for (size_t i = 0; i < vec->size; i++) {
    if (vec->data[i].alto_idx == packed) {
      return &vec->data[i];
    }
  }
  return NULL;
}

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
struct hacoo_tensor *read_init()
{
  return file_init(stdin);
}

/* Read an entry from stdin */
void read_entry(struct hacoo_tensor *t)
{
  file_entry(t, stdin);
}

/* Read a tensor from a tns file */
struct hacoo_tensor *read_tensor_file(FILE *file)
{
  struct hacoo_tensor *t = file_init(file);

  while(!feof(file)) {
    file_entry(t, file);
  }

  return t;
}

/* Merge this with regular function at later time */
struct hacoo_tensor *read_tensor_file_with_base(FILE *file, int zero_base)
{
  struct hacoo_tensor *t = file_init(file);

  while(!feof(file)) {
    file_entry_with_base(t, file, zero_base);
  }

  return t;
}

/* Initialize a tensor from a file */
struct hacoo_tensor *file_init(FILE *file) {

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

struct hacoo_tensor *read_tensor_file_with_base_fast(FILE *file, int zero_base, size_t nnz_estimate)
{
    char line[4096]; // adjust if tensor has many modes
    unsigned int ndims = 0;
    unsigned int *dims = NULL;
    struct hacoo_tensor *t = NULL;

    // --- Read first line to get tensor dimensions ---
    if (!fgets(line, sizeof(line), file)) {
        fprintf(stderr, "Error: Failed to read tensor dimensions.\n");
        return NULL;
    }

    // Count number of dimensions
    char *p = line;
    while (*p) {
        if (*p == ' ') ndims++;
        p++;
    }
    ndims++; // last number

    dims = (unsigned int *)malloc(ndims * sizeof(unsigned int));
    if (!dims) {
        fprintf(stderr, "Error: Failed to allocate dims.\n");
        return NULL;
    }

    // Parse dimensions
    p = line;
    for (unsigned int i = 0; i < ndims; i++) {
        dims[i] = strtoul(p, &p, 10);
    }

    // --- Compute number of buckets ---
    size_t nbuckets = 1ULL << (size_t)ceil(log2((double)nnz_estimate / 0.7));

    // Allocate tensor
    t = hacoo_alloc(ndims, dims, nbuckets, LOAD);
    free(dims);

    // Allocate reusable arrays for indices
    int index_int[ndims];
    unsigned int index_u[ndims];

    // --- Read tensor entries ---
   while (fgets(line, sizeof(line), file)) {
    if (line[0] == '\0' || line[0] == '\n' || line[0] == '#') continue;
    
    char *p = line;
    unsigned int parsed = 0;
    for (unsigned int i = 0; i < ndims; i++) {
        index_int[i] = (int) strtol(p, &p, 10);
        parsed++;
    }

    double value = strtod(p, &p);
    parsed++;

    if (parsed != ndims + 1) {
        fprintf(stderr, "Error: Line does not have expected %u entries: %s\n", ndims + 1, line);
        continue;
    }

    // adjust for 1-based indexing
    if (!zero_base) {
        for (unsigned int i = 0; i < ndims; i++) index_int[i]--;
    }

    for (unsigned int i = 0; i < ndims; i++) index_u[i] = (unsigned int) index_int[i];

    hacoo_set(t, index_u, value);
    }


    return t;
}


/* Read an entry from a file */
void file_entry(struct hacoo_tensor *t, FILE *file) {

  double value;
  unsigned int *index = (unsigned int *) malloc(t->ndims * sizeof(unsigned int));
  if (!index) {
    fprintf(stderr, "Error: Failed to allocate memory for index array.\n");
    return;
  }

  /* read the index */
  for (int i = 0; i < t->ndims; i++) {
    if (feof(file))
      return;
    fscanf(file, "%u", &index[i]);
  }

  /* read the value */
  if (feof(file))
    return;
  fscanf(file, "%lf", &value);

  /* insert the value */
  hacoo_set(t, index, value);
}

/* Merge this with above at later time */
void file_entry_with_base(struct hacoo_tensor *t, FILE *file, int zero_base) {
  
  double value;
  int *index = (int *) malloc(t->ndims * sizeof(int));
  if (!index) {
    fprintf(stderr, "Error: Failed to allocate memory for index array.\n");
    return;
  }

  /* read the index */
  for (int i = 0; i < t->ndims; i++) {
    if (feof(file))
      return;
    fscanf(file, "%u", &index[i]);
  }

  /*if indexes are one-based, like FROSTT tensors, subtract 1*/
  if (!zero_base) {
    
    for (int i = 0; i < t->ndims; i++) {
      if (index[i] == 0) {
        fprintf(stderr, "Error: Tensor uses base-1 indexing but has index 0 in mode %d\n", i);
        free(index);
        exit(EXIT_FAILURE);
      }
      index[i] -= 1;
    }
  }

  unsigned int *uindex = (unsigned int *) malloc(t->ndims * sizeof(unsigned int));
  if (!uindex) {
    fprintf(stderr, "Error: Failed to allocate memory for unsigned index array.\n");
    free(index);
    return;
  }

  for (int i = 0; i < t->ndims; i++) {
    if (index[i] < 0) {
      fprintf(stderr, "Error: Negative index after base adjustment at mode %d: %d\n", i, index[i]);
      free(index);
      free(uindex);
      return;
    }
    uindex[i] = (unsigned int) index[i];
  }

  if (fscanf(file, "%lf", &value) != 1) {
    fprintf(stderr, "Error: Failed to read value\n");
    free(index);
    free(uindex);
    return;
  }

  /* insert the value */
  hacoo_set(t, uindex, value);

  free(index);
  free(uindex);
}

/* Print out information about the tensor */
void print_status(struct hacoo_tensor *t) {
  printf("nnz: %d nbuckets: %d\n", t->nnz, t->nbuckets);
}

void print_tensor(struct hacoo_tensor *t)
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

            printf("0x%llx: ", b->alto_idx);

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
double frobenius_norm(struct hacoo_tensor *t)
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

/*Debugging print functions */
/* Print the nth nonzero element in the tensor */
/*
void print_nth_nonzero(struct hacoo_tensor *t, int n) {
  if (n < 0 || n >= t->nnz) {
      printf("Invalid nonzero index.\n");
      return;
  }

  int count = 0;
  unsigned int index[t->ndims];

  // Iterate over all buckets 
  for (int i = 0; i < t->nbuckets; i++) {
    bucket_vector *vec = &t->buckets[i];
    for (size_t j = 0; j < vec->size; j++) {
      struct hacoo_bucket *b = &vec->data[j];
      if (count == n) {
          hacoo_extract_index(t->buckets[b], t->ndims, index);
          printf("Nonzero %d:\n", n);
          printf("0x%llx: ", t->buckets[b]->morton);
          for (int j = 0; j < t->ndims; j++) {
              printf("%u ", index[j]);
          }
          printf("%f\n", t->buckets[b]->value);
          return;
      }
      count++;
    }
  }

  printf("Nonzero element not found (this should not happen).\n");
}
*/

/* Print the contents of a specific bucket in the tensor */
/*
void print_bucket(struct hacoo_tensor *t, int bucket_index) {
    if (bucket_index < 0 || bucket_index >= t->nbuckets) {
        printf("Invalid bucket index.\n");
        return;
    }

    struct hacoo_bucket *b = t->buckets[bucket_index];
    if (!b) {
        printf("Bucket %d is empty.\n", bucket_index);
        return;
    }

    printf("\nBucket %d:\n=============\n", bucket_index);
    unsigned int index[t->ndims];

    for (; b; b = b->next) {
        hacoo_extract_index(b, t->ndims, index);
        printf("0x%llx: ", b->morton);
        for (int j = 0; j < t->ndims; j++) {
            printf("%u ", index[j]);
        }
        printf("%f\n", b->value);
    }
}
*/

/* Print the contents of a specific bucket */
/*
void print_bucket_from_ptr(struct hacoo_bucket *b, unsigned int ndims) {
    if (!b) {
        printf("Bucket is empty.\n");
        return;
    }

    printf("\nBucket:\n=============\n");
    unsigned int index[ndims];

    for (; b; b = b->next) {
        hacoo_extract_index(b, ndims, index);
        printf("0x%llx: ", b->morton);
        for (unsigned int j = 0; j < ndims; j++) {
            printf("%u ", index[j]);
        }
        printf("%f\n", b->value);
    }
}
*/
