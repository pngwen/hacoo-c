#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Function to construct a 3x3x3 tensor with rank-2
void construct_tensor(double tensor[3][3][3]) {
    // Define factor matrices A, B, C (3x2 each)
    double A[3][2] = {
        {1.0, 0.5},
        {0.8, 0.2},
        {0.3, 0.7}
    };

    double B[3][2] = {
        {0.6, 0.9},
        {0.4, 0.1},
        {0.7, 0.3}
    };

    double C[3][2] = {
        {0.2, 0.8},
        {0.5, 0.6},
        {0.9, 0.4}
    };

    // Compute the tensor as the outer product of A, B, and C
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                tensor[i][j][k] = 0.0;
                for (int r = 0; r < 2; r++) { // Rank-2
                    tensor[i][j][k] += A[i][r] * B[j][r] * C[k][r];
                }
            }
        }
    }
}

// Test the constructed tensor
int main() {
    double tensor[3][3][3];
    construct_tensor(tensor);
    printf("3 3 3\n");

    // Print the tensor
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                printf("%d %d %d %.4f\n", i, j, k, tensor[i][j][k]);
            }
        }
    }

    return 0;
}
