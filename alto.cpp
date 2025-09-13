#include "hacoo.h"
#include "alto.h"
#include "bitops.cpp"
#include <typeinfo>
#include <bitset>
#include <iostream>
#include <cstdint>
#include <cassert>

/*Pack individual index*/
LIT alto_pack_index(const unsigned int *coords,
                    const LIT *ALTO_MASKS,
                    int nmode) {
    LIT alto = 0;
    for (int i = 0; i < nmode; i++) {
        LIT partial = pdep(static_cast<uint64_t>(coords[i]),
                           static_cast<uint64_t>(ALTO_MASKS[i]));
        alto |= partial;

        // Print the original coordinate and its contribution
        /*printf("coords[%d] = %u, ALTO_MASKS[%d] = 0x%llx, partial = 0x%llx\n",
               i,
               coords[i],
               i,
               (unsigned long long) ALTO_MASKS[i],
               (unsigned long long) partial);*/
    }

    // Print the final packed index
    //printf("Packed index = 0x%llx\n", (unsigned long long) alto);

    return alto;
}

/*
LIT alto_pack_index(const unsigned int *coords,
                           const LIT *ALTO_MASKS,
                           int nmode) {
    LIT alto = 0;
    for (int i = 0; i < nmode; i++) {
        alto |= pdep(static_cast<uint64_t>(coords[i]), static_cast<uint64_t>(ALTO_MASKS[i]));
    }
    printf("LTO_MASKS[%d]A = 0x%llx\n", n, ALTO_MASKS[n]);
    return alto;
}*/

// Unpack ALTO indices
void alto_unpack(LIT alto_idx,
                 const LIT* masks,
                 int ndims,
                 unsigned int* out_indices) {
    assert(masks && out_indices);

    //printf("Unpacking index = 0x%llx\n", (unsigned long long) alto_idx);

    for (int m = 0; m < ndims; ++m) {
        LIT extracted = pext(static_cast<uint64_t>(alto_idx),
                             static_cast<uint64_t>(masks[m]));
        out_indices[m] = static_cast<unsigned int>(extracted);

        // Debug print
        /*printf("mask[%d] = 0x%llx, extracted = 0x%llx, out_indices[%d] = %u\n",
               m,
               (unsigned long long) masks[m],
               (unsigned long long) extracted,
               m,
               out_indices[m]);*/
    }
}

//setup packing/linearization with ALTO
// Achieving alto_bits_min requires packing/compression.
void alto_setup(struct hacoo_tensor *at, PackOrder po, ModeOrder mo)
{
    LIT* ALTO_MASKS = (LIT*)calloc(at->ndims, sizeof(LIT));
    assert(ALTO_MASKS);

    int nmode = at->ndims;
    printf("nmode: %d\n", nmode);

    int alto_bits_min = 0, alto_bits_max = 0;
    LIT alto_mask = 0;
    int max_num_bits = 0, min_num_bits = sizeof(IType) * 8;

    MPair* mode_bits = (MPair*)MALLOC(nmode * sizeof(MPair));
    assert(mode_bits);

    // Initial mode values.
    for (int n = 0; n < nmode; ++n) {
        int mbits = (sizeof(IType) * 8) - clz(static_cast<uint64_t>(at->dims[n] - 1));
        mode_bits[n].mode = n;
        mode_bits[n].bits = mbits;
        alto_bits_min += mbits;
        max_num_bits = std::max(max_num_bits, mbits);
        min_num_bits = std::min(min_num_bits, mbits);
        printf("num_bits for mode-%d=%d\n", n+1, mbits);
    }

#ifdef ALT_PEXT
    //Simple prefix sum
    at->mode_pos[0] = 0;
    for (int n = 1; n < nmode; ++n) {
        at->mode_pos[n] = at->mode_pos[n-1] + mode_bits[n-1].bits;
    }
#endif

    printf("Ordinal Type:        int%lu\n", 8*sizeof(IType));
    printf("Sparse Value Type:   %s%lu\n", typeid(ValType) == typeid(double) || typeid(ValType) == typeid(float) ? "FP" : "int", 8*sizeof(ValType));

    alto_bits_max = max_num_bits * nmode;

    //printf("range of mode bits=[%d %d]\n", min_num_bits, max_num_bits);
    printf("alto_bits_min=%d, alto_bits_max=%d\n", alto_bits_min, alto_bits_max);
    assert(alto_bits_min <= ((int)sizeof(LIT) * 8));

    //Assuming we use a power-2 data type for ALTO_idx with a minimum size of a byte
    //int alto_bits = pow(2, (sizeof(int) * 8) - __builtin_clz(alto_bits_min));
    // int alto_bits = (int)0x1 << std::max<int>(3, (sizeof(int) * 8) - __builtin_clz(alto_bits_min));
    // printf("alto_bits=%d\n", alto_bits);

    double alto_storage = 0;
    alto_storage = at->nnz * (sizeof(ValType) + sizeof(LIT));
    printf("Alto format storage:    %g Bytes\n", alto_storage);

    // alto_storage = at->nnz * (sizeof(ValType) + (alto_bits >> 3));
    // printf("Alto-power-2 format storage:    %g Bytes\n", alto_storage);

    alto_storage = at->nnz * (sizeof(ValType) + (alto_bits_min >> 3));
    printf("Alto-opt format storage:    %g Bytes\n", alto_storage);

    {//Dilation & shifting.
        int level = 0, shift = 0, inc = 1;

        //Sort modes, if needed.
        if (mo == SHORT_FIRST)
            std::sort(mode_bits, mode_bits + nmode, [](auto& a, auto& b) { return a.bits < b.bits; });
        else if(mo == LONG_FIRST)
            std::sort(mode_bits, mode_bits + nmode, [](auto& a, auto& b) { return a.bits > b.bits; });

        if (po == MSB_FIRST) {
            shift = alto_bits_min - 1;
            inc = -1;
        }

        bool done;
        do {
            done = true;

            for (int n = 0; n < nmode; ++n) {
                if (level < mode_bits[n].bits) {
                    ALTO_MASKS[mode_bits[n].mode] |= (LIT)0x1 << shift;
                    shift += inc;
                    done = false;
                }
            }
            ++level;
        } while (!done);

        assert(level == (max_num_bits+1));
        assert(po == MSB_FIRST ? (shift == -1) : (shift == alto_bits_min));
    }

    for (int n = 0; n < nmode; ++n) {
        at->mode_masks[n] = ALTO_MASKS[n];
        alto_mask |= ALTO_MASKS[n];
//#ifdef ALTO_DEBUG
        printf("ALTO_MASKS[%d] = 0x%llx\n", n, ALTO_MASKS[n]);
//#endif
    }
    at->alto_mask = alto_mask;
//#ifdef ALTO_DEBUG
    printf("alto_mask = 0x%llx\n", alto_mask);
//#endif

    /*printf("Final ALTO_MASKS (binary):\n");
    for (int n = 0; n < nmode; ++n) {
        printf("Mode %d mask = ", n, ALTO_MASKS[n]);
        for (int b = sizeof(LIT)*8 - 1; b >= 0; --b)
            printf("%d", (ALTO_MASKS[n] >> b) & 1);
        printf("\n");
    } */
    free(mode_bits);
    free(ALTO_MASKS);
}

/* set up for for morton linearization */
void morton_setup(struct hacoo_tensor *at)
{
    at->mode_masks = (LIT*)calloc(at->ndims, sizeof(LIT));
    at->alto_mask = 0;
    
    for(int i=0; i<at->ndims; i++) {
        at->mode_masks[i]=0;
        unsigned int idx=at->dims[i];
        while(idx) {
            at->mode_masks[i] = (at->mode_masks[i] << at->ndims) | 1;
            idx >>= 1;
        }
        alto_mask |= at->mode_masks[i];
    }
}