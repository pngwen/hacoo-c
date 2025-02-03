
/*
 * A simple test program for the hacoo sparse tensor library.
 */
#include "CUnit/Basic.h"
#include "hacoo.h"
#include "matrix.h"
#include "mttkrp.h"
#include <CUnit/CUnit.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Set up MTTKRP CUnit Test Suites */
void CUnit_mttkrp_ser();
void CUnit_mttkrp_par();

/* MTTKRP Tests*/
void mttkrp_test_ser(); // Serial MTTKRP
void mttkrp_test_par(); // Parallel MTTKRP

// Global array to store command-line arguments
char **global_argv;
int global_argc;

int main(int argc, char *argv[]) {

/*
  // read the tensor in
  struct hacoo_tensor *t = read_init();
  while (!feof(stdin)) {
    read_entry(t);
    //print_status(t);
  }

  //Print tensor
  print_tensor(t);
  hacoo_free(t);
*/
  // Store the command-line arguments in the global array
  global_argc = argc;
  global_argv = argv;

  CUnit_mttkrp_ser();
  //CUnit_mttkrp_par();
}

/* Time Serial MTTKRP algorithm.
  Pass test if answers in text files match with this library's mttkrp answers.
  Filenames should be passed through command line.
  global_argv[1] - name of file to read factor matrices from
  global_argv[2] - name of file to read mttkrp results from
  */
void mttkrp_test_ser() {

  CU_ASSERT_PTR_NOT_NULL(global_argv);
  CU_ASSERT_EQUAL(global_argc, 3); // Number of arguments
  const char *factor_matrix_filename = global_argv[1];
  const char *mttkrp_filename = global_argv[2];

  // read factor matrices
  matrix_t **factor_matrices = NULL;

  int matrix_count =
      read_matrices_from_file(factor_matrix_filename, &factor_matrices);
  //printf("matrix count: %d\n", matrix_count);
  if (matrix_count == -1) {
    printf("Error reading matrices.\n");
    return; // Error reading matrices
  }

  // read sptensor data
  struct hacoo_tensor *t = read_init();
  while (!feof(stdin)) {
    read_entry(t);
    //print_status(t);
  }

  // print the tensor
  //print_tensor(t);

  // read mttkrp answers for each mode
  matrix_t **mttkrp_ans = NULL;

  matrix_count = read_matrices_from_file(mttkrp_filename, &mttkrp_ans);

  printf("matrix count: %d\n", matrix_count);
  if (matrix_count == -1) {
    printf("Error reading matrices.\n");
    return; // Error reading matrices
  }

  //Print factor matrices
  //print_matrices(mttkrp_ans, matrix_count);
  
    //calc only mode 3 mttkrp 
      matrix_t **hacoo_mttkrp =
      (matrix_t **)malloc(sizeof(matrix_t *) * t->ndims);
    printf("Mode-3 MTTKRP: \n");
    printf("HaCOO-C Answer: \n");
    hacoo_mttkrp[2] = mttkrp_serial(t, factor_matrices, 2);
    print_matrix(hacoo_mttkrp[2]);
    //printf("MATLAB Answer:\n");
    //print_matrix(mttkrp_ans[2]);
    //CU_ASSERT(are_matrices_equal(mttkrp_ans[2], hacoo_mttkrp[2]));
    free_matrices(hacoo_mttkrp,matrix_count);
/*
  // use this library's mttkrp to calculate answers
  matrix_t **hacoo_mttkrp =
      (matrix_t **)malloc(sizeof(matrix_t *) * t->ndims);
  for (int i = 0; i < matrix_count; i++) {
    printf("\nMode-%d MTTKRP: \n", i + 1);
    hacoo_mttkrp[i] = mttkrp_serial(t, factor_matrices, i);
  }

  // for every mttkrp answer compare answer with this libarary's mttkrp answer
  for (int i = 0; i < matrix_count; i++) {
    printf("HaCOO-C Answer:\n");
    print_matrix(hacoo_mttkrp[i]);
    printf("MATLAB Answer:\n");
    print_matrix(mttkrp_ans[i]);
    CU_ASSERT(are_matrices_equal(mttkrp_ans[i], hacoo_mttkrp[i]));
  }
*/
  // Free factor matrices
  free_matrices(factor_matrices, matrix_count);

  // Free mttkrp answer matrices
  free_matrices(mttkrp_ans, matrix_count);
  //free_matrices(hacoo_mttkrp, matrix_count);

   // free tensor
  hacoo_free(t);
}

/* verify this library's calculation of MTTKRP.
   Pass test if answers in text files match with this library's mttkrp answers.
   Filenames should be passed through command line.
   global_argv[1] - name of file to read factor matrices from
   global_argv[2] - name of file to read mttkrp results from
   */
void mttkrp_test_par() {

  CU_ASSERT_PTR_NOT_NULL(global_argv);
  CU_ASSERT_EQUAL(global_argc, 3); // Number of arguments

  const char *factor_matrix_filename = global_argv[1];
  const char *mttkrp_filename = global_argv[2];

  // read factor matrices
  matrix_t **factor_matrices = NULL;

  int matrix_count =
      read_matrices_from_file(factor_matrix_filename, &factor_matrices);
  //printf("matrix count: %d\n", matrix_count);
  if (matrix_count == -1) {
    printf("Error reading matrices.\n");
    return; // Error reading matrices
  }

  // read sptensor data
  struct hacoo_tensor *t = read_init();
  while (!feof(stdin)) {
    read_entry(t);
    //print_status(t);
  }

  // print the tensor
  //print_tensor(t);

  // read mttkrp answers for each mode
  matrix_t **mttkrp_ans = NULL;

  matrix_count = read_matrices_from_file(mttkrp_filename, &mttkrp_ans);

  //printf("matrix count: %d\n", matrix_count);
  if (matrix_count == -1) {
    printf("Error reading matrices.\n");
    return; // Error reading matrices
  }

  //Print factor matrices
  //print_matrices(mttkrp_ans, matrix_count);

  // use this library's mttkrp to calculate answers
  matrix_t **hacoo_mttkrp =
      (matrix_t **)malloc(sizeof(matrix_t *) * t->ndims);
  for (int i = 0; i < matrix_count; i++) {
    //printf("\nMode-%d MTTKRP: \n", i + 1);
    hacoo_mttkrp[i] = mttkrp(t, factor_matrices, i);
  }

  // for every mttkrp answer compare answer with this libarary's mttkrp answer
  for (int i = 0; i < matrix_count; i++) {
    //printf("HaCOO-C Answer:\n");
    //print_matrix(hacoo_mttkrp[i]);
    printf("HaCOO-C Answer:\n");
    print_matrix(hacoo_mttkrp[i]);
    printf("MATLAB Answer:\n");
    print_matrix(mttkrp_ans[i]);
    CU_ASSERT(are_matrices_equal(mttkrp_ans[i], hacoo_mttkrp[i]));
  }
  // free tensor
  hacoo_free(t);

  // Free factor matrices
  free_matrices(factor_matrices, matrix_count);

  // Free mttkrp answer matrices
  free_matrices(mttkrp_ans, matrix_count);
  free_matrices(hacoo_mttkrp, matrix_count);
}



  /*Serial MTTKRP CUnit Test*/
void CUnit_mttkrp_ser() {
  CU_initialize_registry(); // Initialize the CUnit registry
  CU_pSuite pSuite = CU_add_suite("Serial MTTKRP Test", 0, 0); // Add a test suite
  CU_add_test(pSuite, "Serial MTTKRP", mttkrp_test_ser);
  CU_basic_run_tests();                        // Run the tests
  CU_cleanup_registry(); // Cleanup after running the tests
}

/*Parallel MTTKRP CUnit Test*/
void CUnit_mttkrp_par() {
  CU_initialize_registry(); // Initialize the CUnit registry
  CU_pSuite pSuite = CU_add_suite("Parallel MTTKRP Test", 0, 0); // Add a test suite
  CU_add_test(pSuite, "Parallel MTTKRP", mttkrp_test_par);
  CU_basic_run_tests();                        // Run the tests
  CU_cleanup_registry(); // Cleanup after running the tests
}