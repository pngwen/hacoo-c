/* File: hacoo.c
 * Purpose: Implementation of the hacoo sparse tensor library.
 */
#include "hacoo.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LOAD 70
#define MIN_BUCKETS 128

/* Helper Function Prototypes */
static void hacoo_free_buckets(struct hacoo_tensor *t);
static unsigned long long hacoo_morton(unsigned int n, unsigned int *index);
static size_t hacoo_bucket_index(struct hacoo_tensor *t,
                                 unsigned long long morton);
static void hacoo_compute_params(struct hacoo_tensor *t);
static struct hacoo_bucket *hacoo_bucket_search(struct hacoo_bucket *b,
                                                unsigned long long morton);
static size_t hacoo_max_bits(unsigned int n);

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
  if ((t->nnz / t->nbuckets) > ((float)t->load/100.0))
  {
    hacoo_rehash(t);
  }
}

void hacoo_rehash(struct hacoo_tensor *t)
{
    // Allocate a new tensor with double the number of buckets
    struct hacoo_tensor *new_tensor = hacoo_alloc(t->ndims, t->dims, t->nbuckets * 2, t->load);
    if (!new_tensor)
    {
        fprintf(stderr, "Error: Failed to allocate new tensor during rehash.\n");
        return;
    }

    printf("Rehashing: new number of buckets: %zu\n", new_tensor->nbuckets);

    struct hacoo_bucket *cur;
    unsigned int *index = (unsigned int *)malloc(sizeof(unsigned int) * t->ndims);
    if (!index)
    {
        fprintf(stderr, "Error: Failed to allocate memory for index array during rehash.\n");
        hacoo_free(new_tensor);
        return;
    }

    // Loop over all buckets in the current tensor
    for (size_t i = 0; i < t->nbuckets; i++)
    {
        cur = t->buckets[i];
        while (cur)
        {
            // Extract the index from the current bucket
            hacoo_extract_index(cur, t->ndims, index);

            // Rehash the entry into the new tensor
            hacoo_set(new_tensor, index, cur->value);

            // Move to the next bucket in the chain
            struct hacoo_bucket *next = cur->next;
            free(cur); // Free the current bucket
            cur = next;
        }
    }

    // Free the old tensor's bucket array
    free(index);
    free(t->buckets);

    // Replace the old tensor's fields with the new tensor's fields
    t->buckets = new_tensor->buckets;
    t->nbuckets = new_tensor->nbuckets;
    t->sx = new_tensor->sx;
    t->sy = new_tensor->sy;
    t->sz = new_tensor->sz;

    // Free the temporary tensor structure but not its buckets
    new_tensor->buckets = NULL;
    hacoo_free(new_tensor);
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

/* extract the index from a bucket */
void hacoo_extract_index(struct hacoo_bucket *b, unsigned int n,
                         unsigned int *index)
{
  size_t max_bits = hacoo_max_bits(n);

  // Initialize the index array
  for (unsigned int i = 0; i < n; i++)
  {
    index[i] = 0;
  }
  // De-interleave the morton code bits into the index array
  for (unsigned int bit = 0; bit < max_bits; bit++)
  {
    for (unsigned int i = 0; i < n; i++)
    {
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

static unsigned long long hacoo_morton(unsigned int n, unsigned int *index)
{
  unsigned long long m = 0;
  size_t max_bits = hacoo_max_bits(n);

  for (unsigned int bit = 0; bit < max_bits; bit++)
  {
    for (unsigned int i = 0; i < n; i++)
    {
      m |= ((index[i] >> bit) & 1) << (bit * n + i);
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
  size_t b1 = sizeof(unsigned long long) * 8 / n;
  size_t b2 = sizeof(unsigned int) * 8;

  return b1 < b2 ? b1 : b2;
}

/*Allocate a new hacoo bucket. Its next points to NULL if a next bucket does not
 * exist. */
struct hacoo_bucket *hacoo_new_bucket()
{
  struct hacoo_bucket *b;
  b = calloc(1, sizeof(struct hacoo_bucket));
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
  int index[t->ndims];
  double value;

  /* read the index */
  for (int i = 0; i < t->ndims; i++) {
    if (feof(stdin))
      return;
    fscanf(file, "%u", &index[i]);
    //printf("index: %d\n",index[i]);
  }

  //if FROSTT tensor, subtract 1 from everything
  for(int i=0;i<t->ndims;i++) {
    index[i] = index[i]-1;
  }

  /* read the value */
  if (feof(stdin))
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

/*
int save_hacoo_tensor_to_file(const struct hacoo_tensor *tensor, const char *filename) {
    if (tensor == NULL || filename == NULL) {
        fprintf(stderr, "Error: Invalid tensor or filename.\n");
        return -1;
    }

    FILE *file = fopen(filename, "wb");
    if (file == NULL) {
        perror("Error opening file for writing");
        return -1;
    }

    // Step 1: Write the number of dimensions
    fwrite(&tensor->ndims, sizeof(unsigned int), 1, file);

    // Step 2: Write the dimensions array
    fwrite(tensor->dims, sizeof(unsigned int), tensor->ndims, file);

    // Step 3: Write the number of non-zero entries
    fwrite(&tensor->nnz, sizeof(unsigned int), 1, file);

    size_t *index = (size_t *)malloc(sizeof(size_t) * tensor->ndims);
    // Step 4: Iterate through buckets and save indices and values
    for (int b = 0; b < tensor->nbuckets; b++) {
        struct hacoo_bucket *cur = tensor->buckets[b];
        while (cur != NULL) {

            //TODO: extract indices morton before writing
             hacoo_extract_index(cur, tensor->ndims, index);
            // Write indices (array of size ndims)
            fwrite(index, sizeof(size_t), tensor->ndims, file);

            // Write value (double)
            fwrite(&cur->value, sizeof(double), 1, file);

            cur = cur->next;
        }
    }

    free(index);
    fclose(file);
    return 0; // Success
}

struct hacoo_tensor *load_hacoo_tensor_from_file(const char *filename) {
    if (filename == NULL) {
        fprintf(stderr, "Error: Invalid filename.\n");
        return NULL;
    }

    FILE *file = fopen(filename, "rb");
    if (file == NULL) {
        perror("Error opening file for reading");
        return NULL;
    }

    // Step 1: Read the number of dimensions
    unsigned int ndims;
    if (fread(&ndims, sizeof(unsigned int), 1, file) != 1) {
        fprintf(stderr, "Error reading number of dimensions.\n");
        fclose(file);
        return NULL;
    }

    // Step 2: Read the dimensions array
    unsigned int *dims = (unsigned int *)malloc(sizeof(unsigned int) * ndims);
    if (dims == NULL) {
        fprintf(stderr, "Error allocating memory for dimensions.\n");
        fclose(file);
        return NULL;
    }
    if (fread(dims, sizeof(unsigned int), ndims, file) != ndims) {
        fprintf(stderr, "Error reading dimensions array.\n");
        free(dims);
        fclose(file);
        return NULL;
    }

    // Step 3: Read the number of non-zero entries
    unsigned int nnz;
    if (fread(&nnz, sizeof(unsigned int), 1, file) != 1) {
        fprintf(stderr, "Error reading number of non-zero entries.\n");
        free(dims);
        fclose(file);
        return NULL;
    }

    // Step 4: Allocate memory for the HACOO tensor
    struct hacoo_tensor *tensor = (struct hacoo_tensor *)malloc(sizeof(struct hacoo_tensor));
    if (tensor == NULL) {
        fprintf(stderr, "Error allocating memory for HACOO tensor.\n");
        free(dims);
        fclose(file);
        return NULL;
    }

    tensor->ndims = ndims;
    tensor->dims = dims;
    tensor->nnz = nnz;
    tensor->nbuckets = 1024; // Choose a default bucket size (adjust as needed)
    tensor->buckets = (struct hacoo_bucket **)calloc(tensor->nbuckets, sizeof(struct hacoo_bucket *));
    if (tensor->buckets == NULL) {
        fprintf(stderr, "Error allocating memory for HACOO buckets.\n");
        free(tensor);
        free(dims);
        fclose(file);
        return NULL;
    }

    // Step 5: Read each non-zero entry and add it to the tensor
    for (unsigned int i = 0; i < nnz; i++) {
        // Read the indices
        unsigned int *indices = (unsigned int *)malloc(sizeof(unsigned int) * ndims);
        if (indices == NULL) {
            fprintf(stderr, "Error allocating memory for indices.\n");
            // Free all allocated resources
            for (int b = 0; b < tensor->nbuckets; b++) {
                struct hacoo_bucket *cur = tensor->buckets[b];
                while (cur) {
                    struct hacoo_bucket *next = cur->next;
                    free(cur->idx);
                    free(cur);
                    cur = next;
                }
            }
            free(tensor->buckets);
            free(tensor);
            free(dims);
            fclose(file);
            return NULL;
        }
        if (fread(indices, sizeof(unsigned int), ndims, file) != ndims) {
            fprintf(stderr, "Error reading indices.\n");
            free(indices);
            fclose(file);
            return NULL;
        }

        // Read the value
        double value;
        if (fread(&value, sizeof(double), 1, file) != 1) {
            fprintf(stderr, "Error reading value.\n");
            free(indices);
            fclose(file);
            return NULL;
        }

        // Add the non-zero entry to the tensor (hashing by first index as an example)
        unsigned int bucket = indices[0] % tensor->nbuckets;
        struct hacoo_bucket *new_bucket = (struct hacoo_bucket *)malloc(sizeof(struct hacoo_bucket));
        if (new_bucket == NULL) {
            fprintf(stderr, "Error allocating memory for HACOO bucket.\n");
            free(indices);
            fclose(file);
            return NULL;
        }

        new_bucket->indices = indices;
        new_bucket->value = value;
        new_bucket->next = tensor->buckets[bucket];
        tensor->buckets[bucket] = new_bucket;
    }

    fclose(file);
    return tensor;
}
*/