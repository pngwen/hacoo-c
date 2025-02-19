/* File: hacoo.c
 * Purpose: Implementation of the hacoo sparse tensor library.
 */
#include "hacoo.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define LOAD 70
#define MIN_BUCKETS 128

/* Helper Function Prototypes */
static void hacoo_free_buckets(struct hacoo_tensor *t);
static void free_buckets(struct hacoo_bucket **buckets, int nbuckets);
static uint64_t hacoo_morton(unsigned int n, unsigned int *index);
static size_t hacoo_bucket_index(struct hacoo_tensor *t,
                                 unsigned long long morton);
static void hacoo_compute_params(struct hacoo_tensor *t);
static struct hacoo_bucket *hacoo_bucket_search(struct hacoo_bucket *b,
                                                unsigned long long morton);
static size_t hacoo_max_bits(unsigned int n);

#include <stdio.h>

/* Print the nth nonzero element in the tensor */
void print_nth_nonzero(struct hacoo_tensor *t, int n) {
    if (n < 0 || n >= t->nnz) {
        printf("Invalid nonzero index.\n");
        return;
    }

    int count = 0;
    struct hacoo_bucket *b;
    unsigned int index[t->ndims];

    /* Iterate over all buckets */
    for (int i = 0; i < t->nbuckets; i++) {
        for (b = t->buckets[i]; b; b = b->next) {
            if (count == n) {
                hacoo_extract_index(b, t->ndims, index);
                printf("Nonzero %d:\n", n);
                printf("0x%llx: ", b->morton);
                for (int j = 0; j < t->ndims; j++) {
                    printf("%u ", index[j]);
                }
                printf("%f\n", b->value);
                return;
            }
            count++;
        }
    }

    printf("Nonzero element not found (this should not happen).\n");
}

/* Print the contents of a specific bucket in the tensor */
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

/* Print the contents of a specific bucket */
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

/* Allocation and deallocation functions */
struct hacoo_tensor *hacoo_alloc(unsigned int ndims, unsigned int *dims,
                                 size_t nbuckets, unsigned int load)
{
  struct hacoo_tensor *t = malloc(sizeof(struct hacoo_tensor));

  /* handle allocation error */
  if (t == NULL)
  {
    goto error;
  }

  /* initialize tensor fields */
  t->ndims = ndims;
  t->dims = calloc(1, sizeof(unsigned int) * ndims);
  if (!t->dims)
  {
    goto error;
  }
  memcpy(t->dims, dims, sizeof(unsigned int) * ndims);
  t->nbuckets = nbuckets;
  t->load = load;
  t->buckets = calloc(nbuckets, sizeof(struct hacoo_bucket *));
  if (!t->buckets)
  {
    goto error;
  }
  hacoo_compute_params(t);
  return t;

error:
  if (t)
  {
    hacoo_free(t);
  }
  return NULL;
}

void hacoo_free(struct hacoo_tensor *t)
{
  if (t->dims)
  {
    free(t->dims);
  }
  if (t->buckets)
  {
    hacoo_free_buckets(t);
  }
  free(t);
}

/* Access functions */
void hacoo_set(struct hacoo_tensor *t, unsigned int *index, double value)
{
  //fprintf(stderr, "In hacoo set: recieved index: %u %u %u %u\n",index[0], index[1],index[2],index[3]);
  unsigned long long morton = hacoo_morton(t->ndims, index);
  size_t i = hacoo_bucket_index(t, morton);
  struct hacoo_bucket *b;

  // see if it already exists in the bucket
  b = hacoo_bucket_search(t->buckets[i], morton);

  // If bucket is unoccupied, create new bucket
  if (!b)
  {
    b = hacoo_new_bucket(); // Create a new bucket
    if (t->buckets[i] == NULL)
    {
      // If the bucket at index i is empty, place the new bucket there
      t->buckets[i] = b;
    }
    else
    {
      // If there is already a bucket chain, append to the end of the chain
      struct hacoo_bucket *last = t->buckets[i];
      while (last->next)
      {
        last = last->next;
      }
      last->next = b; // Append the new bucket to the end of the chain
    }
    b->morton = morton; // Set Morton code on the new bucket
    t->nnz++;           // Increment the number of non-zero elements
  }

  b->value = value; // Set the value for the bucket

  /* Check if the load limit has been exceeded*/
  if (t->nbuckets > 0 && ((double)t->nnz / (double)t->nbuckets) > ((double)t->load / 100.0)) {
    hacoo_rehash(&t);
    if (t == NULL) {  // Ensure rehash was successful
      fprintf(stderr, "Rehash failed, exiting.\n");
      return;
    }
}

}

void hacoo_rehash(struct hacoo_tensor **t)
{
  //Allocate dummy tensor with new parameters
  // Step 1: Allocate a new tensor with twice the number of buckets
  struct hacoo_tensor *dummy = hacoo_alloc((*t)->ndims, (*t)->dims, (*t)->nbuckets * 2, (*t)->load);
  hacoo_compute_params(dummy);

  if (dummy== NULL) {
      fprintf(stderr, "Failed to allocate dummy tensor during rehash.\n");
      return;
  }

    // Allocate a new bucket array with double the number of buckets
    struct hacoo_bucket **new_buckets = (struct hacoo_bucket **)calloc((*t)->nbuckets * 2, sizeof(struct hacoo_bucket *));
    if (!new_buckets) {
        fprintf(stderr, "Error: Failed to allocate new buckets during rehash.\n");
        return;
    }

    unsigned int *index = (unsigned int *)malloc(sizeof(unsigned int) * (*t)->ndims);
    if (!index)
    {
        fprintf(stderr, "Error: Failed to allocate memory for index array during rehash.\n");
        free(new_buckets);
        return;
    }

    struct hacoo_bucket *cur;
    struct hacoo_bucket *b;
    int nnz=0; //ensure new nnz count is == current nnz

    // Loop over all buckets in the current tensor
    for (size_t i = 0; i < (*t)->nbuckets; i++) {
      //printf("i: %d\n",i);
      cur = (*t)->buckets[i];
      while (cur) {
        //make new bucket
        b = hacoo_new_bucket(); // Create a new bucket
        b->morton = cur->morton;
        b->value = cur->value;

        // Compute the new bucket index
        size_t j = hacoo_bucket_index(dummy, cur->morton);

        if (new_buckets[j] == NULL) {
          // If the bucket at index j is empty, place the new bucket there
          new_buckets[j] = b;
        }else {
          // If there is already a bucket chain, append to the end of the chain
          struct hacoo_bucket *last = new_buckets[j];
          while (last->next)
          {
            last = last->next;
          }
          last->next = b; // Append the new bucket to the end of the chain
        }
        nnz++;           // Increment the number of non-zero elements
        //printf("inserted nnz %d\n",nnz);
        cur = cur->next;
      }
    }

  //ensure everything copied correctly
  if((*t)->nnz != nnz) {
    printf("Something went wrong. Only %d nnz copied when there were originally %d nnz.\n",nnz,(*t)->nnz);
  }

  // Copy the computed parameters from dummy to t before freeing it
  (*t)->sx = dummy->sx;
  (*t)->sy = dummy->sy;
  (*t)->sz = dummy->sz;

  struct hacoo_bucket **old_buckets = (*t)->buckets; // Save reference to old buckets

  // Update tensor with new buckets and bucket count
  (*t)->buckets = new_buckets;
  (*t)->nbuckets *= 2;

  // Free the old bucket array
  free_buckets(old_buckets,(*t)->nbuckets / 2);
  hacoo_free(dummy);
  free(index);

  //printf("Rehashing completed. tensor has %d nnz and %d buckets.\n",(*t)->nnz,(*t)->nbuckets);
}


double hacoo_get(struct hacoo_tensor *t, unsigned int *index)
{
  unsigned long long morton = hacoo_morton(t->ndims, index);
  unsigned int i = hacoo_bucket_index(t, morton);
  struct hacoo_bucket *b = hacoo_bucket_search(t->buckets[i], morton);
  if (b)
  {
    return b->value;
  }

  return 0.0;
}

/* Extract the index from a bucket */
void hacoo_extract_index(struct hacoo_bucket *b, unsigned int n,
                         unsigned int *index)
{
    size_t max_bits = hacoo_max_bits(n);

    // Initialize the index array
    for (unsigned int i = 0; i < n; i++) {
        index[i] = 0;
    }

    // De-interleave the Morton code bits into the index array
    for (unsigned int bit = 0; bit < max_bits; bit++) {
        for (unsigned int i = 0; i < n; i++) {
            index[i] |= ((b->morton >> (bit * n + i)) & 1) << bit;
        }
    }
}

/* Helper function implementations. */
static void hacoo_free_buckets(struct hacoo_tensor *t)
{
  struct hacoo_bucket *cur, *next;
  int i;

  for (i = 0; i < t->nbuckets; i++)
  {
    cur = t->buckets[i];
    while (cur)
    {
      next = cur->next;
      free(cur);
      cur = next;
    }
  }

  free(t->buckets);
}

/* free only buckets */
static void free_buckets(struct hacoo_bucket **buckets, int nbuckets)
{
  struct hacoo_bucket *cur, *next;
  int i;

  for (i = 0; i < nbuckets; i++)
  {
    cur = buckets[i];
    while (cur)
    {
      next = cur->next;
      free(cur);
      cur = next;
    }
  }

  free(buckets);
}

// Encodes index, returns morton code
// n: ndimensions
static uint64_t hacoo_morton(unsigned int n, unsigned int *index)
{
    uint64_t m = 0;
    size_t max_bits = hacoo_max_bits(n);

    for (unsigned int bit = 0; bit < max_bits; bit++) {
        for (unsigned int i = 0; i < n; i++) {
            m |= ((uint64_t)((index[i] >> bit) & 1)) << (bit * n + i);
        }
    }

    return m;
}

static size_t hacoo_bucket_index(struct hacoo_tensor *t,
                                 unsigned long long morton)
{
  unsigned long long hash = morton;

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

static struct hacoo_bucket *hacoo_bucket_search(struct hacoo_bucket *b,
                                                unsigned long long morton)
{
  for (; b; b = b->next)
  {
    if (b->morton == morton)
    {
      return b;
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

/*Allocate a new hacoo bucket. Its next points to NULL if a next bucket does not
 * exist. */
struct hacoo_bucket *hacoo_new_bucket()
{
  struct hacoo_bucket *b;
  b = calloc(1, sizeof(struct hacoo_bucket));
  if (!b) {
    fprintf(stderr, "Error: Failed to allocate memory for new bucket.\n");
    return NULL;
  }
  b->next = NULL;
  b->morton = 0;
  b->value = 0;

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
  unsigned int *dims = malloc(count * sizeof(unsigned int));
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
  free(dims);

  return t;
}

/* Read an entry from a file */
void file_entry(struct hacoo_tensor *t, FILE *file) {

  double value;
  unsigned int *index = malloc(t->ndims * sizeof(unsigned int));
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

  /*if indexes are base-1, like FROSTT tensors,
  then subtract 1 from everything*/
  /*if(t->base == 1) {
    for (int i = 0; i < t->ndims; i++) {
      index[i] -= index[i];
    }
  }*/

  /* read the value */
  if (feof(file))
    return;
  fscanf(file, "%lf", &value);

  /* insert the value */
  hacoo_set(t, index, value);
}

/* Print out information about the tensor */
void print_status(struct hacoo_tensor *t) {
  printf("nnz: %d nbuckets: %d\n", t->nnz, t->nbuckets);
}

/* Print the tensor hash table with COO listings */
void print_tensor(struct hacoo_tensor *t)
{
  struct hacoo_bucket *b;
  int index[t->ndims];

  for (int i = 0; i < t->nbuckets; i++) {
    // get the bucket pointer
    if (!t->buckets[i])
      continue;

    // print the bucket heading
    printf("\nBucket %d\n=============\n", i);
    for (b = t->buckets[i]; b; b = b->next) {
      hacoo_extract_index(b, t->ndims, index);
      printf("0x%lx: ", b->morton);
      for (int j = 0; j < t->ndims; j++) {
        printf("%u ", index[j]);
      }
      printf("%f\n", b->value);
    }
  }
}
