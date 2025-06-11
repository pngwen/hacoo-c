/* Define matrix_t type- SPLATT does this. Should work nice like a numpy matrix.
Returns a matrix from an array-like object.

Example:
matrix([[1, 2],
        [3, 4]])
*/

#ifndef MATRIX_H
#define MATRIX_H
#include <stdint.h> //for uint32_t which is what SPLATT uses
#include <stdlib.h>

typedef struct matrix {
  unsigned int rows;
  unsigned int cols;
  double **vals;
} matrix_t;

/* Create matrix of all zeros */
matrix_t *new_matrix(unsigned int n_rows, unsigned int n_cols);

/* Generate a random matrix of a given size and value range */
matrix_t* new_random_matrix(size_t rows, size_t cols, double min_value, double max_value);

/* Make a new matrix from existing data */
matrix_t *array_to_matrix(double *data, unsigned int n_rows, unsigned int n_cols);

/* Copy contents of one matrix to a newly allocated matrix */
matrix_t *copy_matrix(matrix_t *m);

/* Copy matrix in place. */
void copy_matrix_to(matrix_t *dest, matrix_t *src);

/* Copy multiple matrices to newly allocated matrices */
matrix_t** copy_matrices(matrix_t **originals, size_t num_matrices);

/* Print values of a matrix */
void print_matrix(matrix_t *m);

void print_matrices(matrix_t **matrices, int num_matrices);

/* free matrix */
matrix_t *free_matrix(matrix_t *m);

/* Read matrix from text file */
matrix_t *matrix_read_init();

/* Read matrix data from a space-delimited file into an array of matrices */
int read_matrices_from_file(const char* filename, matrix_t ***matrices);

// Function to free memory allocated for matrices
void free_matrices(matrix_t **matrices, int matrix_count);

// Check if 2 matrices are equal. returns 1 if considered equal, 0 otherwise.
int are_matrices_equal(const matrix_t *m1, const matrix_t *m2);
        
//check for equality within a specific tolerance
int are_equal(double a, double b);

/* Add two matrices RES = A + B. C must be allocated to the correct dimensions */
void add_matrix(matrix_t *res, matrix_t *a, matrix_t *b);

/* Add two matrices RES = A + B over only a specific column. C must be allocated to the correct dimensions */
void add_matrix_column(matrix_t *res, matrix_t *a, matrix_t *b, int col_idx);

/* Subtract two matrices RES = A - B. C must be allocated to the correct dimensions */
void sub_matrix(matrix_t *res, matrix_t *a, matrix_t *b);

/* Perform the matrix multiplication RES = A * B */
void mul_matrix(matrix_t *res, matrix_t *a, matrix_t *b);

/* Perform the matrix multiplication RES = A' * B */
void mul_transpose_matrix(matrix_t *res, matrix_t *a, matrix_t *b);

/* Multiply each element in the matrix by a scalar */
void scale_matrix(matrix_t *m, double scalar);

/* Calculate the inverse of the matrix RES = A^-1 */
void invert_matrix(matrix_t *res, matrix_t *a);

/* Fill in the identity matrix to an existing matrix */
void fill_identity_matrix(matrix_t *m);

/* Fill a matrix with a number */
void fill_matrix(matrix_t *m, double val);

/* basic test for matrix functions */
void matrix_test();

/* Print 1-D array */
void print_array(void *arr, int size, char type);

void locate_zeroes(double *arr, int size);

void print_matrix_column(matrix_t *matrix, int col_index);

double matrix_frobenius_norm(matrix_t *m);

#endif
