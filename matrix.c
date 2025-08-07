#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>  // For fabs()
#include <string.h>
#include <cblas.h> //for OpenBLAS

//Define acceptable margin of error
#define EPSILON 1.0e-2


matrix_t *new_matrix(unsigned int n_rows, unsigned int n_cols) {
  matrix_t *matrix = (matrix_t *)malloc(sizeof(matrix_t));
  matrix->rows = n_rows;
  matrix->cols = n_cols;
  double **vals = (double **)malloc(sizeof(double *) * n_rows);
  vals[0] = (double *)calloc(n_rows * n_cols, sizeof(double)); // Allocate memory for the first row
  for (int x = 1; x < n_rows; x++) {
    vals[x] = vals[x - 1] + n_cols; // Point subsequent rows to the same memory block
  }
  matrix->data = vals[0];
  matrix->vals = vals;
  return matrix;
}

/* Generate a random matrix of a given size and value range */
matrix_t* new_random_matrix(size_t rows, size_t cols, double min_value, double max_value) {
    matrix_t *random_matrix = new_matrix(rows, cols);
    static int seeded = 0;

    // Seed the random number generator
    if(!seeded) {
        seeded=1;
        srand((unsigned int)time(NULL));
    }

    // Fill the matrix with random values in the specified range [min_value, max_value]
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            random_matrix->vals[i][j] = min_value + (rand() / (double)RAND_MAX) * (max_value - min_value);
        }
    }

    return random_matrix;
}

matrix_t *array_to_matrix(double *data, unsigned int n_rows, unsigned int n_cols) {
  matrix_t *matrix = new_matrix(n_rows, n_cols);
  for (int x = 0; x < n_rows; x++) {
    for (int y = 0; y < n_cols; y++) {
      matrix->vals[x][y] = data[n_cols * x + y];
    }
  }
  return matrix;
}

/* Copy an existing matrix to a newly allocated matrix */
matrix_t* copy_matrix(matrix_t *original) {
    // Create a new matrix with the same dimensions
    matrix_t *copy = new_matrix(original->rows, original->cols);
    memcpy(copy->data, original->data, original->rows * original->cols * sizeof(double));

    return copy;
}

/* Copies the contents of an existing matrix src into a pre-allocated matrix dest */
void copy_matrix_to(matrix_t *dest, matrix_t *src) {
    // Ensure the dest has the same dimensions as the src
    if (src->rows != dest->rows || src->cols != dest->cols) {
        fprintf(stderr, "Error: Dimensions of src and dest do not match.\n");
        return;
    }

    memcpy(dest->data, src->data, src->rows * src->cols * sizeof(double));
}

/* Copy multiple matrices to newly allocated matrices */
matrix_t** copy_matrices(matrix_t **originals, size_t num_matrices) {
    matrix_t **copies = malloc(num_matrices * sizeof(matrix_t *));
    if (!copies) {
        return NULL; // Allocation failed
    }

    for (size_t k = 0; k < num_matrices; k++) {
        matrix_t *original = originals[k];
        if (!original) {
            copies[k] = NULL;
            continue;
        }

        matrix_t *copy = new_matrix(original->rows, original->cols);
        if (!copy) {
            // Clean up previously allocated matrices
            for (size_t j = 0; j < k; j++) {
                free_matrix(copies[j]);
            }
            free(copies);
            return NULL;
        }

        // Use memcpy for the contiguous data array
        memcpy(copy->data, original->data, original->rows * original->cols * sizeof(double));

        copies[k] = copy;
    }

    return copies;
}


void print_matrix(matrix_t *m) {
  //printf("Matrix: (%d x % d)\n", m->rows, m->cols);
  for (int x = 0; x < m->rows; x++) {
    if (x != 0) {
      printf("%s", "\n");
    }
    for (int y = 0; y < m->cols; y++) {
      printf("%f\t", m->vals[x][y]);
    }
  }
  printf("\n");
}

/*print multiple matrices */
void print_matrices(matrix_t **matrices, int num_matrices) {
  for (int i = 0; i < num_matrices; i++) {
    printf("Matrix %d (%d x %d):\n", i + 1, matrices[i]->rows,
           matrices[i]->cols); // print with dimensions
    matrix_t *m = matrices[i]; // Get the matrix pointer

    // Print each matrix
    for (int x = 0; x < m->rows; x++) {
      if (x != 0) {
        printf("%s", "\n");
      }
      for (int y = 0; y < m->cols; y++) {
        printf("%f\t", m->vals[x][y]);
      }
    }
    printf("\n\n"); // Adding a space between matrices for clarity
  }
}

/* Free a single matrix */
void free_matrix(matrix_t *m) {
    if (!m) { return; }
    if(m->data) free(m->data);
    if(m->vals) free(m->vals);
    free(m);
}

// Function to read matrices from a file into an array of matrix_t pointers
int read_matrices_from_file(const char *filename, matrix_t ***matrices) {
  FILE *file = fopen(filename, "r");
  if (!file) {
    perror("Cannot open file");
    return -1;
  }

  int matrix_count = 0;
  matrix_t **matrix_array = NULL; // Pointer to array of matrix_t pointers

  while (!feof(file)) {
    int rows, cols;

    // Read dimensions of the matrix
    if (fscanf(file, "%d %d", &rows, &cols) != 2) {
      break; // End of file or invalid data
    }

    // Allocate memory for the new array of matrix_t pointers
    matrix_t **new_matrix_array =
        (matrix_t **)malloc((matrix_count + 1) * sizeof(matrix_t *));
    if (!new_matrix_array) {
      perror("Memory allocation failed");
      fclose(file);
      return -1;
    }

    // Copy the previous matrices into the new array
    if (matrix_count > 0) {
      for (int i = 0; i < matrix_count; i++) {
        new_matrix_array[i] = matrix_array[i];
      }
      free(matrix_array); // Free the old array after copying
    }

    // Allocate memory for the new matrix
    new_matrix_array[matrix_count] = new_matrix(rows, cols);
    if (!new_matrix_array[matrix_count]) {
      perror("Memory allocation failed");
      fclose(file);
      free(new_matrix_array); // Free newly allocated memory
      return -1;
    }


    // Read the matrix data
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        if (fscanf(file, "%lf", &new_matrix_array[matrix_count]->vals[i][j]) != 1) {
          fprintf(stderr, "Error reading matrix data at row %d, col %d (matrix %d)\n", i, j, matrix_count);
          fclose(file);
          free(new_matrix_array);
          return -1;
        }
      }
    }

    // Update matrix_count and matrix_array
    matrix_count++;
    matrix_array = new_matrix_array;
  }

  fclose(file);
  *matrices = matrix_array;
  return matrix_count; // Return the number of matrices read
}

// Free an array of matrices
void free_matrices(matrix_t **matrices, int matrix_count) {
  for (int i = 0; i < matrix_count; i++) {
      free_matrix(matrices[i]); // Free each matrix
  }
  free(matrices);
}

// Compare two matrices for equality within a given tolerance
int are_matrices_equal(const matrix_t *m1, const matrix_t *m2) {
  if (m1->rows != m2->rows || m1->cols != m2->cols) {
    return 0; // return false
  }

  for (int i = 0; i < m1->rows; i++) {
    for (int j = 0; j < m1->cols; j++) {
      if (!are_equal(m1->vals[i][j], m2->vals[i][j])) {
        return 0;
      }
    }
  }

  return 1; // return true
}

//check for equality within a specific tolerance
int are_equal(double a, double b) {
    return fabs(a - b) < EPSILON * fmax(fabs(a), fabs(b));
}
/* Test to read matrix from text file */
int read_matrix_test(const char *filename) {

  matrix_t **matrices = NULL;

  // Read matrices from the file
  int matrix_count = read_matrices_from_file(filename, &matrices);
  printf("matrix count: %d\n", matrix_count);
  if (matrix_count == -1) {
    return 1; // Error reading matrices
  }

  // Print the matrices and their dimensions
  for (int i = 0; i < matrix_count; i++) {
    printf("Matrix %d (%d x %d):\n", i + 1, matrices[i]->rows,
           matrices[i]->cols);
    print_matrix(matrices[i]);
    printf("%s\n", "----------------------");
  }

  // Free memory
  free_matrices(matrices, matrix_count);

  return 0;
}

// for testing
void matrix_test() {

  double a[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  matrix_t *m1 = array_to_matrix(a, 4, 3);
  print_matrix(m1);
  free_matrix(m1);
}

void add_matrix(matrix_t *res, matrix_t *a, matrix_t *b) {
  for (int i = 0; i < a->rows; i++) {
    for (int j = 0; j < a->cols; j++) {
      res->vals[i][j] = a->vals[i][j] + b->vals[i][j];
    }
  }
}

//Add values only over a specifix column
void add_matrix_column(matrix_t *res, matrix_t *a, matrix_t *b, int col_idx) {
  // Ensure the column index is valid
  if (col_idx < 0 || col_idx >= a->cols) {
    printf("Invalid column index\n");
    return;
  }

  for (int i = 0; i < a->rows; i++) {
    // Add only over the specified column
    res->vals[i][col_idx] = a->vals[i][col_idx] + b->vals[i][col_idx];
  }
}


void sub_matrix(matrix_t *res, matrix_t *a, matrix_t *b) {
  for (int i = 0; i < a->rows; i++) {
    for (int j = 0; j < a->cols; j++) {
      res->vals[i][j] = a->vals[i][j] - b->vals[i][j];
    }
  }
}

void mul_matrix(matrix_t *res, matrix_t *a, matrix_t *b) {
    if (a->cols != b->rows || res->rows != a->rows || res->cols != b->cols) {
        fprintf(stderr, "Matrix dimension mismatch in mul_matrix\n");
        return;
    }

    cblas_dgemm(
        CblasRowMajor,      // Row-major storage
        CblasNoTrans,       // No transpose on A
        CblasNoTrans,       // No transpose on B
        a->rows,            // M
        b->cols,            // N
        a->cols,            // K
        1.0,                // Alpha
        a->data,            // A
        a->cols,            // lda
        b->data,            // B
        b->cols,            // ldb
        0.0,                // Beta
        res->data,          // C
        res->cols           // ldc
    );
}

/*Old mul matrix*/
/*
void mul_matrix(matrix_t *res, matrix_t *a, matrix_t *b)
{
    matrix_t *tmp = new_matrix(a->rows, b->cols);

    for (int i = 0; i < a->rows; i++)
    {
        for (int j = 0; j < b->cols; j++)
        {
            tmp->vals[i][j] = 0.0;
            for (int k = 0; k < a->cols; k++)
            {
                tmp->vals[i][j] += a->vals[i][k] * b->vals[k][j];
            }
        }
    }

    // Copy the result to res
    copy_matrix_to(res, tmp);
    free_matrix(tmp);
}
*/

void mul_transpose_matrix(matrix_t *res, matrix_t *a, matrix_t *b)
{
    matrix_t *tmp = new_matrix(a->cols, b->cols);

    for (int i = 0; i < a->cols; i++)
    {
        for (int j = 0; j < b->cols; j++)
        {
            tmp->vals[i][j] = 0.0;
            for (int k = 0; k < a->rows; k++)
            {
                tmp->vals[i][j] += a->vals[k][i] * b->vals[k][j];
            }
        }
    }

    // Copy the result to res
    copy_matrix_to(res, tmp);
    free_matrix(tmp);
}


/* Multiply each element in the matrix by a scalar */
void scale_matrix(matrix_t *m, double scalar)
{
    for (int i = 0; i < m->rows; i++)
    {
        for (int j = 0; j < m->cols; j++)
        {
            m->vals[i][j] *= scalar;
        }
    }
}


void invert_matrix(matrix_t *res, matrix_t *a)
{
    if (a->rows != a->cols)
        return;

    int n = a->rows;
    matrix_t *augmented = new_matrix(n, 2 * n);

    // Copy a into left half, identity into right half
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            augmented->vals[i][j] = a->vals[i][j];
        for (int j = n; j < 2 * n; j++)
            augmented->vals[i][j] = (i == j - n) ? 1.0 : 0.0;
    }

    // Gaussian elimination with partial pivoting
    for (int i = 0; i < n; i++) {
        // Find pivot row
        int max_row = i;
        double max_val = fabs(augmented->vals[i][i]);
        for (int k = i + 1; k < n; k++) {
            if (fabs(augmented->vals[k][i]) > max_val) {
                max_val = fabs(augmented->vals[k][i]);
                max_row = k;
            }
        }
        if (max_val < 1e-12) { // Singular matrix
            goto done;
        }
        // Swap rows if needed
        if (max_row != i) {
            double *tmp = augmented->vals[i];
            augmented->vals[i] = augmented->vals[max_row];
            augmented->vals[max_row] = tmp;
        }
        // Normalize pivot row
        double pivot = augmented->vals[i][i];
        for (int j = 0; j < 2 * n; j++)
            augmented->vals[i][j] /= pivot;
        // Eliminate other rows
        for (int k = 0; k < n; k++) {
            if (k == i) continue;
            double factor = augmented->vals[k][i];
            for (int j = 0; j < 2 * n; j++)
                augmented->vals[k][j] -= factor * augmented->vals[i][j];
        }
    }

    // Copy right half to result
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            res->vals[i][j] = augmented->vals[i][j + n];

done:
    free_matrix(augmented);
}

/* Print 1-D Array */
void print_array(void *arr, int size, char type) {
    printf("[");
    for (int i = 0; i < size; i++) {
        switch (type) {
            case 'f':  // double
                printf("%.2f", ((double *)arr)[i]);
                break;
            case 'd':  // int or unsigned int
                printf("%d", ((int *)arr)[i]);
                break;
            case 'u':  // explicitly unsigned int
                printf("%u", ((unsigned int *)arr)[i]);
                break;
            default:
                printf("?");
                break;
        }

        if (i < size - 1) {
            printf(", ");
        }
    }
    printf("]\n");
}

/*Print a speciic matrix column */
void print_matrix_column(matrix_t *matrix, int col_index) {
    if (col_index < 0 || col_index >= matrix->cols) {
        printf("Invalid column index!\n");
        return;
    }

    printf("Column %d: [", col_index);
    for (int i = 0; i < matrix->rows; i++) {
        printf("%.2f", matrix->vals[i][col_index]);
        if (i < matrix->rows - 1) {
            printf(", ");
        }
    }
    printf("]\n");
}


/* Fill in the identity matrix to an existing matrix */
void fill_identity_matrix(matrix_t *m)
{
    for(int i=0; i<m->rows; i++)
    {
        for(int j=0; j<m->cols; j++)
        {
            if(i == j)
            {
                m->vals[i][j] = 1;
            }
            else
            {
                m->vals[i][j] = 0;
            }
        }
    }
}

/* Fill a matrix with a number */
void fill_matrix(matrix_t *m, double val)
{
    for(int i=0; i<m->rows; i++)
    {
        for(int j=0; j<m->cols; j++)
        {
            m->vals[i][j] = val;
        }
    }
}


double matrix_frobenius_norm(matrix_t *m) {
    double sum = 0.0;
    for (int i = 0; i < m->rows; i++)
        for (int j = 0; j < m->cols; j++)
            sum += m->vals[i][j] * m->vals[i][j];
    return sqrt(sum);
}
