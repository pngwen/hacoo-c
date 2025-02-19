#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h> // for time
#include <math.h>  // For fabs()

//Define acceptable margin of error
#define EPSILON 1.0e-2


matrix_t *new_matrix(unsigned int n_rows, unsigned int n_cols) {
  matrix_t *matrix = (matrix_t *)malloc(sizeof(matrix_t));
  matrix->rows = n_rows;
  matrix->cols = n_cols;
  double **vals = (double **)malloc(sizeof(double *) * n_rows);
  for (int x = 0; x < n_rows; x++) {
    vals[x] = (double *)calloc(n_cols, sizeof(double));
  }
  matrix->vals = vals;
  return matrix;
}

/* Generate a random matrix of a given size and value range */
matrix_t* new_random_matrix(size_t rows, size_t cols, double min_value, double max_value) {
    matrix_t *random_matrix = new_matrix(rows, cols);

    // Seed the random number generator
    srand((unsigned int)time(NULL));

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

    // Copy the data from the original matrix to the new one
    for (size_t i = 0; i < original->rows; i++) {
        for (size_t j = 0; j < original->cols; j++) {
            copy->vals[i][j] = original->vals[i][j];
        }
    }

    return copy;
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
matrix_t *free_matrix(matrix_t *m) {
  for (int x = 0; x < m->rows; x++) {
    free(m->vals[x]);
  }
  free(m->vals);
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
    new_matrix_array[matrix_count] = (matrix_t *)malloc(sizeof(matrix_t));
    if (!new_matrix_array[matrix_count]) {
      perror("Memory allocation failed");
      fclose(file);
      free(new_matrix_array); // Free newly allocated memory
      return -1;
    }

    // Set the matrix dimensions
    new_matrix_array[matrix_count]->rows = rows;
    new_matrix_array[matrix_count]->cols = cols;
    new_matrix_array[matrix_count]->vals =
        (double **)malloc(rows * sizeof(double *));
    if (!new_matrix_array[matrix_count]->vals) {
      perror("Memory allocation failed");
      fclose(file);
      free(new_matrix_array); // Free memory on error
      return -1;
    }

    // Allocate memory for each row in the matrix
    for (int i = 0; i < rows; i++) {
      new_matrix_array[matrix_count]->vals[i] =
          (double *)malloc(cols * sizeof(double));
      if (!new_matrix_array[matrix_count]->vals[i]) {
        perror("Memory allocation failed");
        fclose(file);
        free(new_matrix_array); // Free memory on error
        return -1;
      }
    }

    // Read the matrix data
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        if (fscanf(file, "%lf", &new_matrix_array[matrix_count]->vals[i][j]) !=
            1) {
          perror("Error reading matrix data");
          fclose(file);
          free(new_matrix_array); // Free memory on error
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
    for (int j = 0; j < matrices[i]->rows; j++) {
      free(matrices[i]->vals[j]);
    }
    free(matrices[i]->vals);
    free(matrices[i]);
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

void sub_matrix(matrix_t *res, matrix_t *a, matrix_t *b) {
  for (int i = 0; i < a->rows; i++) {
    for (int j = 0; j < a->cols; j++) {
      res->vals[i][j] = a->vals[i][j] - b->vals[i][j];
    }
  }
}

void mul_matrix(matrix_t *res, matrix_t *a, matrix_t *b)
{
    for (int i = 0; i < a->rows; i++)
    {
        for (int j = 0; j < b->cols; j++)
        {
            res->vals[i][j] = 0;
            for (int k = 0; k < a->cols; k++)
            {
                res->vals[i][j] += a->vals[i][k] * b->vals[k][j];
            }
        }
    }
}

void mul_transpose_matrix(matrix_t *res, matrix_t *a, matrix_t *b)
{
    for (int i = 0; i < a->cols; i++)
    {
        for (int j = 0; j < b->cols; j++)
        {
            res->vals[i][j] = 0;
            for (int k = 0; k < a->rows; k++)
            {
                res->vals[i][j] += a->vals[k][i] * b->vals[k][j];
            }
        }
    }
}

void invert_matrix(matrix_t *res, matrix_t *a)
{
    // Check if the matrix is square
    if (a->rows != a->cols)
    {
        //TODO: better error handling
        return;
    }

    // Create a temporary matrix to store the augmented matrix
    matrix_t *augmented = new_matrix(a->rows, 2 * a->cols);

    // Copy the original matrix into the left half of the augmented matrix
    for (int i = 0; i < a->rows; i++)
    {
        for (int j = 0; j < a->cols; j++)
        {
            augmented->vals[i][j] = a->vals[i][j];
        }
    }

    // Fill the right half of the augmented matrix with the identity matrix
    for (int i = 0; i < a->rows; i++)
    {
        for (int j = a->cols; j < 2 * a->cols; j++)
        {
            augmented->vals[i][j] = (i == j - a->cols) ? 1 : 0;
        }
    }

    // Perform row operations to transform the left half of the augmented matrix into the identity matrix
    for (int i = 0; i < a->rows; i++)
    {
        double pivot = augmented->vals[i][i];
        if (pivot == 0)
        {
            // TODO: better error handling
            free_matrix(augmented);
            return;
        }

        // Divide the row by the pivot to make the diagonal element 1
        for (int j = 0; j < 2 * a->cols; j++)
        {
            augmented->vals[i][j] /= pivot;
        }

        // Subtract multiples of the row from the other rows to make the rest of the column zero
        for (int k = 0; k < a->rows; k++)
        {
            if (k == i)
            {
                continue;
            }

            double factor = augmented->vals[k][i];
            for (int j = 0; j < 2 * a->cols; j++)
            {
                augmented->vals[k][j] -= factor * augmented->vals[i][j];
            }
        }
    }

    // Copy the right half of the augmented matrix into the result matrix
    for (int i = 0; i < a->rows; i++)
    {
        for (int j = 0; j < a->cols; j++)
        {
            res->vals[i][j] = augmented->vals[i][j + a->cols];
        }
    }
}

/* Print 1-D Array */

void print_array(double *arr, int size) {
    printf("Array: [");
    for (int i = 0; i < size; i++) {
        printf("%d", arr[i]);
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