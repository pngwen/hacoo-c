//#include "CUnit/Basic.h"
#include "hacoo.h"
#include "matrix.h"
#include "mttkrp.h"
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static unsigned long long hacoo_morton(unsigned int n, unsigned int *index);
/* extract the index from a bucket */
void hacoo_extract_index(struct hacoo_bucket *b, unsigned int n,
                         unsigned int *index);
static size_t hacoo_max_bits(unsigned int n);

int main() {

    const int NDIMS = 4;

    //first few lines of uber tensor
    int indexes[10][4] = {
        {0, 0, 535, 664},
        {0, 0, 538, 688},
        {0, 0, 538, 898},
        {0, 0, 541, 894},
        {0, 0, 541, 897},
        {0, 0, 542, 897},
        {0, 0, 547, 718},
        {0, 0, 553, 734},
        {0, 0, 561, 807},
        {0, 0, 569, 717}
    };

    //to store morton encodings
    unsigned long long morton[10];

    //encode each
    for (int i = 0; i < 10; i++) {
        morton[i] = hacoo_morton(NDIMS, indexes[i]);
        printf("%llu\n",morton[i]);
    }

    return 0;
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

// Encodes index, returns morton code
// n: ndimensions
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

static size_t hacoo_max_bits(unsigned int n)
{
  size_t b1 = sizeof(unsigned long long) * 8 / n;
  size_t b2 = sizeof(unsigned int) * 8;

  return b1 < b2 ? b1 : b2;
}