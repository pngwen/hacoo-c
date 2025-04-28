//#include "CUnit/Basic.h"
#include "hacoo.h"
#include "matrix.h"
#include "mttkrp.h"
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Define function pointer type for MTTKRP */
typedef matrix_t *(*mttkrp_func_t)(struct hacoo_tensor *, matrix_t **, unsigned int);

/* Basic read tensor test */
void read_and_print(int argc, char *argv[]);

/* Set up MTTKRP CUnit Test Suite */
void CUnit_mttkrp();
void verify_mttkrp(int mode);
void mttkrp_test_function(void);

/* Calculate mttkrp over all modes */
matrix_t **get_mttkrp_results(struct hacoo_tensor *t, matrix_t **factor_matrices, 
  int matrix_count, mttkrp_func_t f, int mode);

// Global array to store command-line arguments
char **global_argv;
int global_argc;

struct hacoo_tensor *global_tensor = NULL;
matrix_t **global_factor_matrices = NULL;
matrix_t **global_mttkrp_answers = NULL;
int global_matrix_count = 0;

// Global variable to store function pointer
mttkrp_func_t selected_mttkrp_func;

// Global current mode for testing
int current_mode = 0;

int main(int argc, char *argv[]) {
   
  // Store command-line arguments for CUnit access
  global_argc = argc;
  global_argv = argv;

  CUnit_mttkrp();

  return 0;
}

/* Verify MTTKRP algorithm.
   Pass test if MATLAB answers match with this library's MTTKRP answers.
*/
void verify_mttkrp(int mode) {
  CU_ASSERT_PTR_NOT_NULL(global_tensor);
  CU_ASSERT_PTR_NOT_NULL(global_factor_matrices);
  CU_ASSERT_PTR_NOT_NULL(global_mttkrp_answers);

  matrix_t *hacoo_mttkrp = selected_mttkrp_func(global_tensor, global_factor_matrices, mode);

  if (are_matrices_equal(global_mttkrp_answers[mode], hacoo_mttkrp)) {
      CU_PASS("MTTKRP over mode succeeded.\n");
  } else {
      CU_FAIL("MTTKRP over mode failed");
      printf("Failure over Mode %d.\n", mode);
      printf("HaCOO-C Answer:\n");
      print_matrix(hacoo_mttkrp);
      printf("MATLAB Answer:\n");
      print_matrix(global_mttkrp_answers[mode]);
  }

  free_matrix(hacoo_mttkrp);
}

/* Test function that runs verify_mttkrp using current_mode */
void mttkrp_test_function(void) {
    verify_mttkrp(current_mode);
}

/* Calculate MTTKRP over all modes & store in array of matrices */
matrix_t **get_mttkrp_results(struct hacoo_tensor *t, matrix_t **factor_matrices, int matrix_count, mttkrp_func_t f, int mode) {
    matrix_t **hacoo_mttkrp = (matrix_t **)malloc(sizeof(matrix_t *) * t->ndims);
    for (int i = 0; i < matrix_count; i++) {
        hacoo_mttkrp[i] = f(t, factor_matrices, i);
    }
    return hacoo_mttkrp;
}

void CUnit_mttkrp() {
  CU_initialize_registry();
  CU_pSuite pSuite = CU_add_suite("MTTKRP Test", 0, 0);

  const int mode = atoi(global_argv[4]);

  if (mode == 0) {
      selected_mttkrp_func = mttkrp_serial;
      printf("Running Serial MTTKRP Test\n");
  } else if (mode == 1) {
      selected_mttkrp_func = mttkrp;
      printf("Running Parallel MTTKRP Test\n");
  } else {
      printf("Invalid mode. Quitting.\n");
      return;
  }

  // Read tensor
  FILE *file = fopen(global_argv[1], "r");
  if (!file) {
      perror("Error opening tensor file");
      return;
  }
  global_tensor = read_tensor_file(file);
  fclose(file);

  // Read factor matrices
  global_matrix_count = read_matrices_from_file(global_argv[2], &global_factor_matrices);
  if (global_matrix_count == -1) {
      printf("Error reading factor matrices.\n");
      return;
  }

  // Read expected MTTKRP answers
  global_matrix_count = read_matrices_from_file(global_argv[3], &global_mttkrp_answers);
  if (global_matrix_count == -1) {
      printf("Error reading MTTKRP answers.\n");
      return;
  }

  int ndims = global_tensor->ndims;

  // Dynamically add a test for each mode
  for (int i = 0; i < ndims; i++) {
      char test_name[50];
      snprintf(test_name, sizeof(test_name), "MTTKRP Mode %d", i);

      current_mode = i; // set the mode
      CU_add_test(pSuite, test_name, mttkrp_test_function);
  }

  CU_basic_run_tests();
  CU_cleanup_registry();

  // Free global memory
  free_matrices(global_factor_matrices, global_matrix_count);
  free_matrices(global_mttkrp_answers, global_matrix_count);
  hacoo_free(global_tensor);
}

/* Read tensor and print it */
void read_and_print(int argc, char *argv[]) {
  FILE *file = fopen(argv[1], "r");

  if (!file) {
    perror("Error opening file");
    return;
  }

  // Read the tensor
  struct hacoo_tensor *t = read_tensor_file(file);
  fclose(file);

  // Print tensor
  print_tensor(t);

  // Free tensor
  hacoo_free(t);
}