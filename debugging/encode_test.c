#include "CUnit/Basic.h"
#include "hacoo.h"
#include "matrix.h"
#include "mttkrp.h"
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>  // Include for uint64_t

static uint64_t hacoo_morton(unsigned int n, unsigned int *index);
/* extract the index from a bucket */
void hacoo_extract_index(struct hacoo_bucket *b, unsigned int n,
                         unsigned int *index);
static size_t hacoo_max_bits(unsigned int n);

int main() {

    const int NDIMS = 4;

    // First few lines of uber tensor
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

    // To store Morton encodings
    uint64_t morton[10];

    // Encode each
    for (int i = 0; i < 10; i++) {
        morton[i] = hacoo_morton(NDIMS, indexes[i]);
        printf("%llu\n", (unsigned long long)morton[i]);  // Cast for printing
    }

    return 0;
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

// Encodes index, returns Morton code
// n: number of dimensions
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

static size_t hacoo_max_bits(unsigned int n)
{
    size_t b1 = sizeof(uint64_t) * 8 / n;
    size_t b2 = sizeof(unsigned int) * 8;

    return b1 < b2 ? b1 : b2;
}
