#include <stdio.h>
#include <stdint.h>

uint64_t encode(const int *vals, int n) {
    uint64_t result = 0;
    uint64_t bit = 1;
    int max_bits = 32; // Assuming 32-bit integers
    
    for (int b = 0; b < max_bits; b++) {
        for (int i = 0; i < n; i++) {
            if ((vals[i] >> b) & 1) {
                result |= bit;
            }
            bit <<= 1;
        }
    }
    return result;
}

void decode(uint64_t m, int n, int *l) {
    for (int i = 0; i < n; i++) {
        l[i] = 0;
    }
    
    uint64_t bit = 1;
    for (int b = 0; b < 32; b++) {
        for (int i = 0; i < n; i++) {
            if (m & bit) {
                l[i] |= (1 << b);
            }
            bit <<= 1;
        }
    }
}

int main() {
    int indexes[][4] = {
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
    
    int n = 4;
    for (int i = 0; i < 10; i++) {
        uint64_t m = encode(indexes[i], n);
        printf("Morton code: %llu\n", (long long)m);
        
        int decoded[4];
        decode(m, n, decoded);
        printf("Decoded: [");
        for (int j = 0; j < n; j++) {
            printf("%d%s", decoded[j], j == n - 1 ? "" : ", ");
        }
        printf("]\n");
    }
    return 0;
}
