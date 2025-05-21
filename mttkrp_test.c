//#include "CUnit/Basic.h"
#include "hacoo.h"
#include "matrix.h"
#include "mttkrp.h"
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void read_and_print(int argc, char *argv[]);

/* Function pointer type for MTTKRP */
typedef matrix_t *(*mttkrp_func_t)(struct hacoo_tensor *, matrix_t **, unsigned int);

/* Command-line arguments */
char **global_argv;
int global_argc;

/* Function pointer to MTTKRP function */
mttkrp_func_t selected_mttkrp_func;

/* Global data loaded once for tests */
struct hacoo_tensor *global_tensor = NULL;
matrix_t **global_factors = NULL;
matrix_t **global_mttkrp_expected = NULL;
int global_matrix_count = 0;

/*CUnit initialization & cleanup*/
int suite_cleanup(void);
int suite_init(void);

/* CUnit test to verify if computed answers equal matlab's*/
void verify_mttkrp();
void CUnit_verify_mttkrp(); 
matrix_t **get_mttkrp_results(struct hacoo_tensor *t, matrix_t **factor_matrices, int matrix_count, mttkrp_func_t f);

/*Compare algorithm speeds*/
void CUnit_mttkrp_algorithm_comp();

/* Main function */
int main(int argc, char *argv[]) {
    // Force immediate output
    setvbuf(stdout, NULL, _IONBF, 0); // disable stdout buffering
    setvbuf(stderr, NULL, _IONBF, 0); // disable stderr buffering
    global_argc = argc;
    global_argv = argv;
    
    //CUnit_verify_mttkrp();
    CUnit_mttkrp_algorithm_comp();

    return 0;
}

//run MTTKRP algorithm speed comparison
void CUnit_mttkrp_algorithm_comp() {
    CU_initialize_registry();
    CU_pSuite pSuite = CU_add_suite("MTTKRP Modes", 0, 0);
    suite_init();

    const int alg = atoi(global_argv[4]);
    if (alg == 0) {
        selected_mttkrp_func = mttkrp_serial;
        printf("Running Serial MTTKRP Test\n");
    } else if (alg == 1) {
        selected_mttkrp_func = mttkrp;
        printf("Running Parallel MTTKRP Test\n");
    } else {
        printf("Invalid algorithm option. Quitting.\n");
        CU_cleanup_registry();
        return;
    }

    double total_time = 0.0;

    for (int i = 0; i < global_tensor->ndims; i++) {
        struct timespec start, end;
        clock_gettime(CLOCK_MONOTONIC, &start);

        matrix_t *computed = selected_mttkrp_func(global_tensor, global_factors, i);

        clock_gettime(CLOCK_MONOTONIC, &end);

        double duration = (end.tv_sec - start.tv_sec) + 
                          (end.tv_nsec - start.tv_nsec) / 1e9;
        total_time += duration;

        printf("Mode %d MTTKRP Time: %.9f seconds\n", i, duration);

        free_matrix(computed);
    }

    double avg_time = total_time / global_tensor->ndims;
    printf("Average MTTKRP Time across %d modes: %.9f seconds\n", global_tensor->ndims, avg_time);

    CU_basic_run_tests();
    suite_cleanup();
    CU_cleanup_registry();
}

//check this library's mttkrp with matlab's answer over 1 mode
void verify_mttkrp_mode(int mode) {
    CU_ASSERT_PTR_NOT_NULL(global_tensor);
    CU_ASSERT_PTR_NOT_NULL(global_factors);
    CU_ASSERT_PTR_NOT_NULL(global_mttkrp_expected);
    CU_ASSERT(global_matrix_count > 0);

    matrix_t *computed = selected_mttkrp_func(global_tensor, global_factors, mode);
    matrix_t *expected = global_mttkrp_expected[mode];

    if (are_matrices_equal(expected, computed)) {
        CU_PASS("MTTKRP mode test passed.");
    } else {
        CU_FAIL("MTTKRP mode test failed.");
        printf("Failure in Mode %d.\n", mode);
        //printf("HaCOO-C Answer:\n");
        //print_matrix(computed);
        //printf("MATLAB Answer:\n");
        //print_matrix(expected);
    }

    free_matrix(computed);
}

/* Verify if computed answers are equal to matlab's */
void CUnit_verify_mttkrp() {
    CU_initialize_registry();
    suite_init();
    CU_pSuite pSuite = CU_add_suite("MTTKRP Test", 0, 0);

    // Set the MTTKRP implementation
    const int alg = atoi(global_argv[4]);
    if (alg == 0) {
        selected_mttkrp_func = mttkrp_serial;
        printf("Running Serial MTTKRP Test\n");
    } else if (alg== 1) {
        selected_mttkrp_func = mttkrp;
        printf("Running Parallel MTTKRP Test\n");
    } else {
        printf("Invalid algorithm option. Quitting.\n");
        CU_cleanup_registry();
        return;
    }

    CU_add_test(pSuite, "MTTKRP", verify_mttkrp);
    CU_basic_run_tests();
    suite_cleanup();
    CU_cleanup_registry();
}

/* Verify MTTKRP computation with MATLAB answers*/
void verify_mttkrp() {
    CU_ASSERT_PTR_NOT_NULL(global_tensor);
    CU_ASSERT_PTR_NOT_NULL(global_factors);
    CU_ASSERT_PTR_NOT_NULL(global_mttkrp_expected);
    CU_ASSERT(global_matrix_count > 0);

    matrix_t **computed = get_mttkrp_results(global_tensor, global_factors, global_matrix_count, selected_mttkrp_func);

    for (int i = 0; i < global_matrix_count; i++) {
        if (are_matrices_equal(global_mttkrp_expected[i], computed[i])) {
            CU_PASS("MTTKRP over mode succeeded.");
        } else {
            CU_FAIL("MTTKRP over mode failed.");
            printf("Failure over Mode %d.\n", i + 1);
            printf("HaCOO-C Answer:\n");
            print_matrix(computed[i]);
            printf("MATLAB Answer:\n");
            print_matrix(global_mttkrp_expected[i]);
        }
    }

    free_matrices(computed, global_matrix_count);
}

/* Compute MTTKRP over all modes */
matrix_t **get_mttkrp_results(struct hacoo_tensor *t, matrix_t **factor_matrices, int matrix_count, mttkrp_func_t f) {
    matrix_t **results = (matrix_t **)malloc(sizeof(matrix_t *) * t->ndims);
    for (int i = 0; i < matrix_count; i++) {
        results[i] = f(t, factor_matrices, i);
    }
    return results;
}

/* Suite initialization: read all input files */
int suite_init(void) {
    const char *tensor_filename = global_argv[1];
    const char *factor_filename = global_argv[2];
    const char *mttkrp_filename = global_argv[3];

    // Read tensor
    FILE *file = fopen(tensor_filename, "r");
    if (!file) {
        perror("Error opening tensor file");
        exit(1);
    }
    global_tensor = read_tensor_file(file);
    fclose(file);
    if (!global_tensor) return 1;

    // Read factor matrices
    global_matrix_count = read_matrices_from_file(factor_filename, &global_factors);
    if (global_matrix_count == -1) {
        fprintf(stderr, "Error reading factor matrices.\n");
        exit(1);
    }

    // Read expected MTTKRP results
    int expected_count = read_matrices_from_file(mttkrp_filename, &global_mttkrp_expected);
    if (expected_count == -1 || expected_count != global_matrix_count) {
        fprintf(stderr, "Error reading expected MTTKRP results.\n");
        return 1;
    }

    return 0;
}

/* Suite cleanup: free all loaded data */
int suite_cleanup(void) {
    if (global_tensor) {
        hacoo_free(global_tensor);
        global_tensor = NULL;
    }
    if (global_factors) {
        free_matrices(global_factors, global_matrix_count);
        global_factors = NULL;
    }
    if (global_mttkrp_expected) {
        free_matrices(global_mttkrp_expected, global_matrix_count);
        global_mttkrp_expected = NULL;
    }
    global_matrix_count = 0;

    return 0;
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