/* Implementation  of MTTKRP via Sparse Tensor-Vector products (according to
 * Algorithm 1 in the SPLATT paper) */

#include "hacoo.h"
#include "matrix.h"

/* Perform MTTKRP on sparse HaCOO tensor t */
matrix_t *mttkrp(struct hacoo_tensor *t, matrix_t **u, unsigned int n);

/* Serial version of MTTKRP */
matrix_t *mttkrp_serial(struct hacoo_tensor *h, matrix_t **u, unsigned int n);

void mttkrp_test(struct hacoo_tensor *t);

void resizeIntArray(unsigned int** arr, int originalSize);

void resizeDoubleArray(double** arr, int originalSize);