/* HaCOO mmplementation  of MTTKRP via Sparse Tensor-Vector products (according to
 * Algorithm 1 in the SPLATT paper) */

#include "hacoo.h"
#include "matrix.h"

matrix_t *hacoo_mttkrp_debug(struct hacoo_tensor *h, matrix_t **u, unsigned int n);

/* Perform MTTKRP on sparse HaCOO tensor t */
matrix_t *hacoo_mttkrp(struct hacoo_tensor *t, matrix_t **u, unsigned int n);

/* Serial version of MTTKRP */
matrix_t *hacoo_mttkrp_serial(struct hacoo_tensor *h, matrix_t **u, unsigned int n);

void hacoo_mttkrp_test(struct hacoo_tensor *t);