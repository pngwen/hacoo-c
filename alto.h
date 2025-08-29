/* Encoding functions from ALTO */

#ifndef ALTO_H
#define ALTO_H

#include "hacoo.h"

//ALTO currently supports multiple options:
// 1) packing (LSB first or MSB first)
// 2) mode order within a group of bits (natural, shortest first, longest first)
typedef enum PackOrder_ { LSB_FIRST, MSB_FIRST } PackOrder;
typedef enum ModeOrder_ { SHORT_FIRST, LONG_FIRST, NATURAL } ModeOrder;

struct MPair {
    int mode;
    int bits;
};

/* Set up linearization scheme*/
void alto_setup(struct hacoo_tensor *at, PackOrder po, ModeOrder mo);

/* Pack index */
LIT alto_pack_index(const unsigned int *coords, const LIT *ALTO_MASKS, int nmode);

/* Unpack index */
void alto_unpack(unsigned int alto_idx, const unsigned int* masks, int ndims, unsigned int* out_indices);

#endif