#include "hacoo.h"
#include "alto.h"
#include "bitops.cpp"
#include <typeinfo>
#include <bitset>
#include <iostream>
#include <cstdint>

/*Pack individual index*/
LIT alto_pack_index(const unsigned int *coords,
                           const LIT *ALTO_MASKS,
                           int nmode) {
    LIT alto = 0;
    for (int i = 0; i < nmode; i++) {
        alto |= pdep(static_cast<uint64_t>(coords[i]), static_cast<uint64_t>(ALTO_MASKS[i]));
    }
    return alto;
}

//this needs to be changed to 64
#include <cstdint>
#include <vector>
#include <cassert>

// `ALTO_MASKS` contains one mask per mode
// `nmode` is the number of modes
// Returns a vector of indices, one per mode
/*std::vector<uint64_t> alto_unpack(uint64_t alto_idx, const LIT* ALTO_MASKS, int nmode) {
    std::vector<uint64_t> indices(nmode, 0);

    for (int m = 0; m < nmode; ++m) {
        // Extract the bits corresponding to this mode
        indices[m] = pext(alto_idx, ALTO_MASKS[m]);
    }

    return indices;
}*/

// Unpack ALTO indices
#include <immintrin.h> // for _pext_u64
#include <cassert>
#include <cstdint>

void alto_unpack(unsigned int alto_idx, const unsigned int* masks, int ndims, unsigned int* out_indices) {
    assert(masks && out_indices);
    for (int m = 0; m < ndims; ++m) {
        out_indices[m] = static_cast<unsigned int>(pext(static_cast<uint64_t>(alto_idx), masks[m]));
    }
}



/* bulk pack indexes
void
alto_pack(struct hacoo_tensor *at, int nprtn)
{
    /*
    //uint64_t ticks;
    double wtime_s, wtime;

    assert(spt->nmodes <= MAX_NUM_MODES);
    int nmode = spt->nmodes;
    IType nnz = spt->nnz;

    struct hacoo_tensor * _at = (AltoTensor<LIT>*)MALLOC(sizeof(AltoTensor<LIT>));
    assert(_at);

    _at->nmode = nmode;
    _at->nprtn = nprtn;
    _at->nnz = nnz;
    _at->is_oidx = false;

    _at->dims = (IType*)MALLOC(nmode * sizeof(IType));
    assert(_at->dims);
    memcpy(_at->dims, spt->dims, nmode * sizeof(IType));

    _at->mode_masks = (LIT*)MALLOC(nmode * sizeof(LIT));
    assert(_at->mode_masks);
    memset(_at->mode_masks, 0, nmode * sizeof(LIT));

#ifdef ALT_PEXT
    _at->mode_pos = (int*)MALLOC(nmode * sizeof(int));
    assert(_at->mode_pos);
#endif

    _at->idx = (LIT*)MALLOC(nnz * sizeof(LIT));
    assert(_at->idx);

    _at->vals = (ValType*)MALLOC(nnz * sizeof(ValType));
    assert(_at->vals);

    _at->prtn_ptr = (IType*)MALLOC((nprtn + 1) * sizeof(IType));
    assert(_at->prtn_ptr);

    _at->prtn_intervals = (Interval*)MALLOC(static_cast<IType>(nprtn) * static_cast<IType>(nmode) * sizeof(Interval));
    assert(_at->prtn_intervals);

    _at->cr_masks = (LIT*)MALLOC(nmode * sizeof(LIT));
    assert(_at->cr_masks);
    memset(_at->cr_masks, 0, nmode * sizeof(LIT));

#ifdef OPT_ALTO
    _at->prtn_id = (LPType*)MALLOC(nprtn * sizeof(LPType));
    assert(_at->prtn_id);
    memset(_at->prtn_id, 0, nprtn * sizeof(LPType));

    _at->prtn_mask = (LPType*)MALLOC(nprtn * sizeof(LPType));
    assert(_at->prtn_mask);
    memset(_at->prtn_mask, 0, nprtn * sizeof(LPType));

    _at->prtn_mode_masks = (LPType*)MALLOC(nprtn * nmode * sizeof(LPType));
    assert(_at->prtn_mode_masks);
    memset(_at->prtn_mode_masks, 0, nprtn * nmode * sizeof(LPType));
#endif


    //Setup the linearization scheme.
    //ticks = ReadTSC();
    wtime_s = omp_get_wtime();
    alto_setup(_at, LSB_FIRST, SHORT_FIRST);
    //wtime = ElapsedTime (ReadTSC() - ticks);
    wtime = omp_get_wtime() - wtime_s;
    printf("ALTO: setup time = %f (s)\n", wtime);

    //local buffer
    LIT ALTO_MASKS[MAX_NUM_MODES];
    for (int n = 0; n < nmode; ++n) {
        ALTO_MASKS[n] = _at->mode_masks[n];
    }

    //Linearization
    wtime_s = omp_get_wtime();
    #pragma omp parallel for
    for (IType i = 0; i < nnz; i++) {
        LIT alto = 0;

        _at->vals[i] = (ValType)spt->vals[i];
        for (int j = 0; j < nmode; j++) {
            alto |= pdep(spt->cidx[j][i], ALTO_MASKS[j]);
        }
        _at->idx[i] = alto;
#ifdef TEST_ALTO
        for (int j = 0; j < nmode; j++) {
            IType mode_idx = 0;
            mode_idx = pext(alto, ALTO_MASKS[j]);
            assert(mode_idx == spt->cidx[j][i]);
        }
#endif
    }
    wtime = omp_get_wtime() - wtime_s;
    printf("ALTO: Linearization time = %f (s)\n", wtime);

    //Sort the nonzeros based on their line position.
    wtime_s = omp_get_wtime();
    sort_alto(_at);
    wtime = omp_get_wtime() - wtime_s;
    printf("ALTO: sort time = %f (s)\n", wtime);

#ifdef ALT_PEXT
    //Re-encode the ALTO index.

    //local buffer
    int ALTO_POS[MAX_NUM_MODES];
    for (int n = 0; n < nmode; ++n) {
        ALTO_POS[n] = _at->mode_pos[n];
    }

    wtime_s = omp_get_wtime();
    #pragma omp parallel for schedule(static)
    for (IType i = 0; i < nnz; i++) {
        LIT index = _at->idx[i];
        LIT new_index = 0;
        for (int n = 0; n < nmode; ++n) {
            LIT mode_idx = pext(index, ALTO_MASKS[n]);
            new_index |= (mode_idx << ALTO_POS[n]);
        }
        _at->idx[i] = new_index;
    }
    //Update the mode masks to match num_bits.
    for (int n = 0; n < nmode; ++n) {
        int num_bits = (sizeof(IType) * 8) - clz(_at->dims[n] - 1);
        _at->mode_masks[n] = ((1 << num_bits) - 1);
    }
    wtime = omp_get_wtime() - wtime_s;
    printf("ALTO: reorder time = %f (s)\n", wtime);
#ifdef ALTO_DEBUG
    for (int n = 0; n < nmode; n++) {
        printf("ALTO_MASKS[%d] = 0x%llx, pos=%d\n", n, _at->mode_masks[n], _at->mode_pos[n]);
    }
#endif
#endif

    wtime_s = omp_get_wtime();
    //Generate oidx data if needed.
    for (int n = 0; n < nmode; ++n) {
        IType fib_reuse = _at->nnz / _at->dims[n];
        printf("ALTO: fib_reuse[%d]=%llu\n", n, fib_reuse);
        if (fib_reuse <= MIN_FIBER_REUSE) {
            _at->is_oidx = true;
            break;
        }
    }//modes

    if (_at->is_oidx) {
        printf("Tensor requires output-oriented traversal.\n");
        _at->oidx = (IType**) MALLOC(nmode * sizeof(IType*));
        assert(_at->oidx);
        double oidx_storage = 0;
        for (int n = 0; n < nmode; ++n) {
            IType fib_reuse = _at->nnz / _at->dims[n];
            if (fib_reuse <= MIN_FIBER_REUSE) {
                _at->oidx[n] = (IType*) MALLOC(nnz * sizeof(IType));
                assert(_at->oidx[n]);
                oidx_storage += nnz * sizeof(IType);
                #pragma omp parallel for
                for (IType i = 0; i < nnz; ++i) {
                    _at->oidx[n][i] = i;
                }
 #ifndef ALT_PEXT
                std::sort(_at->oidx[n], (_at->oidx[n]) + nnz, [&](IType x, IType y) {
                    return (pext(_at->idx[x], ALTO_MASKS[n]) < pext(_at->idx[y], ALTO_MASKS[n]));});
#else
                std::sort(_at->oidx[n], (_at->oidx[n]) + nnz, [&](IType x, IType y) {
                    return (pext(_at->idx[x], ALTO_MASKS[n], _at->mode_pos[n]) < pext(_at->idx[y], ALTO_MASKS[n], _at->mode_pos[n]));});
#endif

            } else _at->oidx[n] = NULL;
        }//modes
        printf("Auxiliary format storage:    %g Bytes\n", oidx_storage);
    }
    wtime = omp_get_wtime() - wtime_s;
    printf("ALTO: oidx time = %f (s)\n", wtime);

    //Workload partitioning
    wtime_s = omp_get_wtime();
    prtn_alto(_at, nprtn);
    wtime = omp_get_wtime() - wtime_s;
    printf("ALTO: prtn time = %f (s)\n", wtime);

    *at = _at;
}*/

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

    // Initial mode values
    for (int n = 0; n < nmode; ++n) {
        int mbits = (sizeof(IType) * 8) - clz(static_cast<uint64_t>(at->dims[n] - 1));
        mode_bits[n].mode = n;
        mode_bits[n].bits = mbits;
        alto_bits_min += mbits;
        max_num_bits = std::max(max_num_bits, mbits);
        min_num_bits = std::min(min_num_bits, mbits);
        printf("num_bits for mode-%d=%d\n", n, mbits);
    }

#ifdef ALT_PEXT
    at->mode_pos[0] = 0;
    for (int n = 1; n < nmode; ++n) {
        at->mode_pos[n] = at->mode_pos[n-1] + mode_bits[n-1].bits;
    }
#endif

    printf("Ordinal Type:        int%lu\n", 8*sizeof(IType));
    printf("Sparse Value Type:   %s%lu\n",
           typeid(ValType) == typeid(double) || typeid(ValType) == typeid(float) ? "FP" : "int",
           8*sizeof(ValType));

    alto_bits_max = max_num_bits * nmode;
    printf("alto_bits_min=%d, alto_bits_max=%d\n", alto_bits_min, alto_bits_max);
    assert(alto_bits_min <= ((int)sizeof(LIT) * 8));

    // Dilation & shifting
    int level = 0, shift = 0, inc = 1;
    if (mo == SHORT_FIRST)
        std::sort(mode_bits, mode_bits + nmode, [](auto& a, auto& b) { return a.bits < b.bits; });
    else if (mo == LONG_FIRST)
        std::sort(mode_bits, mode_bits + nmode, [](auto& a, auto& b) { return a.bits > b.bits; });

    if (po == MSB_FIRST) {
        shift = alto_bits_min - 1;
        inc = -1;
    }

    bool done;
    do {
        done = true;
        printf("level=%d, shift=%d\n", level, shift);
        for (int n = 0; n < nmode; ++n) {
            if (level < mode_bits[n].bits) {
                ALTO_MASKS[mode_bits[n].mode] |= (LIT)0x1 << shift;

                // Print updated mask in binary
                printf("Updated ALTO_MASKS[%d] = ", n, ALTO_MASKS[n]);
                for (int b = sizeof(LIT)*8 - 1; b >= 0; --b) {
                    printf("%d", (ALTO_MASKS[n] >> b) & 1);
                }
                printf(" (level=%d, shift=%d)\n", level, shift);

                shift += inc;
                done = false;
            }
        }
        ++level;
    } while (!done);


    printf("Final ALTO_MASKS (binary):\n");
    for (int n = 0; n < nmode; ++n) {
        printf("Mode %d mask = ", n, ALTO_MASKS[n]);
        for (int b = sizeof(LIT)*8 - 1; b >= 0; --b)
            printf("%d", (ALTO_MASKS[n] >> b) & 1);
        printf("\n");
    }

    // Copy masks to tensor and compute combined mask
    for (int n = 0; n < nmode; ++n) {
        at->mode_masks[n] = ALTO_MASKS[n];
        alto_mask |= ALTO_MASKS[n];
    }
    at->alto_mask = alto_mask;

    printf("Combined alto_mask = 0x%llx\n", alto_mask);

    free(mode_bits);
    free(ALTO_MASKS);
}
