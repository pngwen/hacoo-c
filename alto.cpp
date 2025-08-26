
#include "hacoo.h"
#include "alto.h"
#include "bitops.hpp"
#include <typeinfo>

/*void
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
void
alto_setup(struct hacoo_tensor *at, PackOrder po, ModeOrder mo)
{
    //LIT ALTO_MASKS[at->ndims] = {}; //initialized to zeros by default
    LIT* ALTO_MASKS = (LIT*)calloc(at->ndims, sizeof(LIT));
    assert(ALTO_MASKS);


    int nmode = at->ndims;
    printf("nmode: %d\n",at->ndims);
    int alto_bits_min = 0, alto_bits_max = 0;
    LIT alto_mask = 0;
    int max_num_bits = 0, min_num_bits = sizeof(IType) * 8;

    MPair* mode_bits = (MPair*)MALLOC(nmode * sizeof(MPair));
    assert(mode_bits);

    //Initial mode values.
    for (int n = 0; n < nmode; ++n) {
        int mbits = (sizeof(IType) * 8) - clz(at->dims[n] - 1);
        mode_bits[n].mode = n;
        mode_bits[n].bits = mbits;
        alto_bits_min += mbits;
        max_num_bits = std::max(max_num_bits, mbits);
        min_num_bits = std::min(min_num_bits, mbits);
        printf("num_bits for mode-%d=%d\n", n + 1, mbits);
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

            printf("here\n");

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
                    printf("ALTO_MASKS[%d]=0x%llx, bits=%d\n", n, ALTO_MASKS[n], mode_bits[n].bits);
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
#ifdef ALTO_DEBUG
        printf("ALTO_MASKS[%d] = 0x%llx\n", n, ALTO_MASKS[n]);
#endif
    }
    at->alto_mask = alto_mask;
#ifdef ALTO_DEBUG
    printf("alto_mask = 0x%llx\n", alto_mask);
#endif
    free(mode_bits);
    free(ALTO_MASKS);
}