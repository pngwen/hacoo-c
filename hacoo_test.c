#include "CUnit/Basic.h"
#include "hacoo.h"
#include "matrix.h"
#include "mttkrp.h"
#include <CUnit/CUnit.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Define function pointer type for MTTKRP */
typedef matrix_t *(*mttkrp_func_t)(struct hacoo_tensor *, matrix_t **, int);

/* Basic read tensor test */
void read_and_print(int argc, char *argv[]);

/* Set up MTTKRP CUnit Test Suite */
void CUnit_mttkrp();
void verify_mttkrp();

/* Calculate mttkrp over all modes */
matrix_t **get_mttkrp_results(struct hacoo_tensor *t, matrix_t **factor_matrices, int matrix_count, mttkrp_func_t f);

// Global array to store command-line arguments
char **global_argv;
int global_argc;

// Global variable to store function pointer
mttkrp_func_t selected_mttkrp_func;

int main(int argc, char *argv[]) {

  //read_and_print(argc,argv);

  // Store command-line arguments for CUnit access
  global_argc = argc;
  global_argv = argv;

  CUnit_mttkrp();

  return 0;
}

/* Verify MTTKRP algorithm.
  Pass test if MATLAB answers match with this library's MTTKRP answers.
*/
void verify_mttkrp() {
    CU_ASSERT_PTR_NOT_NULL(global_argv);
    CU_ASSERT_EQUAL(global_argc, 5); // Expecting 4 arguments + program name

    const char *tensor_filename = global_argv[1];
    const char *factor_matrix_filename = global_argv[2];
    const char *mttkrp_filename = global_argv[3];

    // Read factor matrices
    matrix_t **factor_matrices = NULL;
    int matrix_count = read_matrices_from_file(factor_matrix_filename, &factor_matrices);
    
    if (matrix_count == -1) {
        printf("Error reading factor matrices.\n");
        return;
    }

    FILE *file = fopen(tensor_filename, "r");
    if (!file) {
        perror("Error opening tensor file");
        return;
    }

    // Read the tensor
    struct hacoo_tensor *t = read_tensor_file(file);
    fclose(file);

    // Read expected MTTKRP answers
    matrix_t **mttkrp_ans = NULL;
    matrix_count = read_matrices_from_file(mttkrp_filename, &mttkrp_ans);
    
    if (matrix_count == -1) {
        printf("Error reading MTTKRP matrices.\n");
        return;
    }

    // Compute MTTKRP results using the selected function
    matrix_t **hacoo_mttkrp = get_mttkrp_results(t, factor_matrices, matrix_count, selected_mttkrp_func);

    // Compare computed results with expected answers
    for (int i = 0; i < matrix_count; i++) {
        printf("HaCOO-C Answer:\n");
        print_matrix(hacoo_mttkrp[i]);
        printf("MATLAB Answer:\n");
        print_matrix(mttkrp_ans[i]);
        CU_ASSERT(are_matrices_equal(mttkrp_ans[i], hacoo_mttkrp[i]));
    }

    // Free allocated memory
    free_matrices(factor_matrices, matrix_count);
    free_matrices(mttkrp_ans, matrix_count);
    free_matrices(hacoo_mttkrp, matrix_count);
    hacoo_free(t);
}

/* Calculate MTTKRP over all modes & store in array of matrices */
matrix_t **get_mttkrp_results(struct hacoo_tensor *t, matrix_t **factor_matrices, int matrix_count, mttkrp_func_t f) {
    matrix_t **hacoo_mttkrp = (matrix_t **)malloc(sizeof(matrix_t *) * t->ndims);
    for (int i = 0; i < matrix_count; i++) {
        hacoo_mttkrp[i] = f(t, factor_matrices, i);
    }
    return hacoo_mttkrp;
}

void CUnit_mttkrp() {
    CU_initialize_registry();
    CU_pSuite pSuite = CU_add_suite("MTTKRP Test", 0, 0);

    // Set the function pointer before running the test
    //0 for serial, 1 for parallel
    const int mode = atoi(global_argv[4]);

    if (mode == 0) {
      selected_mttkrp_func = mttkrp_serial;
      printf("Running Serial MTTKRP Test\n");
    } else if (mode == 1) {
      selected_mttkrp_func = mttkrp;
      printf("Running Parallel MTTKRP Test\n");
    } else {
      printf("Invalid mode. Quitting. \n");
      return NULL;
    }

    CU_add_test(pSuite, "MTTKRP", verify_mttkrp);
    CU_basic_run_tests();
    CU_cleanup_registry();
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