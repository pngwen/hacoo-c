#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <getopt.h>
#include <omp.h>
#include <cblas.h>
#include "hacoo.h"
#include "matrix.h"
#include "mttkrp.h"
#include "common.cpp"
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>

/* Function pointer type for MTTKRP */
typedef matrix_t *(*mttkrp_func_t)(struct hacoo_tensor *, matrix_t **, unsigned int);

/* Function pointer to MTTKRP function */
mttkrp_func_t selected_mttkrp_func;

/* Function declarations */
void print_usage(const char *progname);

/* Functions for benchmarking MTTRKP */
int suite_bench_init(const char *tensor_filename, int zero_base, int rank, int nnz);
int generate_factor_matrices();
void CUnit_mttkrp_bench(const char *tensor_file, int alg, int zero_base,
                        int target_mode, int rank, int num_threads, int num_iterations, int nnz);
int suite_cleanup(void);

/* Globals */
struct hacoo_tensor *global_tensor = NULL;
matrix_t **global_factors = NULL;
matrix_t **global_mttkrp_expected = NULL;
int global_matrix_count = 0;
char *global_factor_file = NULL;
char *global_mttkrp_expected_file = NULL;


/* CUnit test to verify if this libary's MTTKRP answers are correct */
int suite_verify_init(const char *tensor_filename, const char *factor_filename, const char *mttkrp_filename, int zero_base, int nnz);
void verify_mttkrp();
void CUnit_verify_mttkrp(const char *tensor_file, const char *factor_file, const char *mttkrp_file, int alg, int zero_base,int nnz);
matrix_t **get_mttkrp_results(struct hacoo_tensor *t, matrix_t **factor_matrices, int matrix_count, mttkrp_func_t f);

void print_usage(const char *progname) {
    printf("Usage: %s [OPTIONS]\n", progname);
    printf("Options:\n");
    printf("  -i or --input          Input tensor file (.tns)\n");
    printf("  -v or --nnz           NUmber of nonzeros\n");
    printf("  -f or --factors        Path to factor matrices\n");
    printf("  -e or --expected       Path to expected MTTKRP answers\n");
    printf("  -z or --zero-based     Assume input tensor is zero-based (default: one-based)\n");
    printf("  -r or --rank           Rank (default: 16)\n");
    printf("  -m or --target-mode    Target mode of tensor (default: all modes)\n");
    printf("  -a or --algorithm      (-2: sequential, default; -1: OpenMP parallel)\n");
    printf("  -b or --bench          Run benchmark mode\n");
    printf("  -d or --dims           Dimensions (I,J,K)\n");
    printf("  -h or --help           Display this help message\n");
    printf("OpenMP options:\n");
    printf("  -t or --number-threads Number of threads (default: 1)      \n");
    printf("\n");
}

/* Main function */
int main(int argc, char *argv[]) {
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

    char *tensor_file = NULL;
    char *dims_str = NULL;
    int algorithm = -2; //default sequential
    int rank = 16;
    int zero_base = 0;
    int target_mode = -1; //default all modes
    int run_bench = 0;
    int num_threads = 1;
	int num_iterations = 1;
	int nnz = 128;

    int opt;
	const char* const short_opt = "hi:v:za:r:m:d:bt:f:e:n:";
	static struct option long_options[] = {
		{"help",        no_argument,       0, 'h'},
		{"input",       required_argument, 0, 'i'},
		{"nnz",         required_argument, 0, 'v'},
		{"factors",     required_argument, 0, 'f'},
		{"expected-mttkrp", required_argument, 0, 'e'},
		{"zero-based",  no_argument,       0, 'z'},
		{"algorithm",   required_argument, 0, 'a'},
		{"rank",        required_argument, 0, 'r'},
		{"target-mode", required_argument, 0, 'm'},
		{"dims",        required_argument, 0, 'd'},
		{"bench",       no_argument,       0, 'b'},
		{"number-threads", required_argument, 0, 't'},
		{"iterations",  required_argument, 0, 'n'},  // NEW
		{0, 0, 0, 0}
	};


    while ((opt = getopt_long(argc, argv, short_opt, long_options, NULL)) != -1) {
        switch (opt) {
            case 'h':
                print_usage(argv[0]);
                exit(0);
            case 'i':
                tensor_file = optarg;
                break;
            case 'v':
                nnz = atoi(optarg);
                break;
            case 'f':
                global_factor_file = optarg;
                break;
            case 'e':
                global_mttkrp_expected_file = optarg;
                break;
            case 'z':
                zero_base = 1;
                break;
			case 'n':
				num_iterations = atoi(optarg);
				if (num_iterations <= 0) {
					fprintf(stderr, "Invalid number of iterations: %s\n", optarg);
					exit(1);
				}
				break;
            case 'r':
                rank = atoi(optarg);
                break;
            case 'm':
                target_mode = atoi(optarg);
                break;
            case 'a':
                algorithm = atoi(optarg);
                break;
            case 'd':
                dims_str = optarg;
                break;
            case 'b':
                run_bench = 1;
                break;
            case 't':
                num_threads = atoi(optarg);
                if (num_threads <= 0) {
                    fprintf(stderr, "Invalid number of threads: %s\n", optarg);
                    exit(1);
                }
                break;
            default:
                print_usage(argv[0]);
                exit(1);
        }
    }

    if (!tensor_file) {
        fprintf(stderr, "Missing required tensor file... exiting\n");
        print_usage(argv[0]);
        exit(1);
    }

    omp_set_num_threads(num_threads);
    openblas_set_num_threads(num_threads);
	if (run_bench) {
		CUnit_mttkrp_bench(tensor_file, algorithm, zero_base, target_mode, rank, num_threads,num_iterations, nnz);
	} else {
		rank = 4;
		printf("Verifying MTTKRP answers.\nTensor: %s\nRank automatically set to 4.\n", tensor_file);
		CUnit_verify_mttkrp(tensor_file, global_factor_file, global_mttkrp_expected_file, algorithm, zero_base,nnz);
	}

    return 0;
}

void CUnit_mttkrp_bench(const char *tensor_file, int alg, int zero_base,
                        int target_mode, int rank, int num_threads, int num_iterations, int nnz) {
    // Initialize CUnit
    CU_initialize_registry();
    
    if (suite_bench_init(tensor_file, zero_base, rank,nnz)) {
        fprintf(stderr, "Suite initialization failed.\n");
        CU_cleanup_registry();
        return;
    }

    // Select MTTKRP implementation
    if (alg == -1) {
        selected_mttkrp_func = mttkrp;
        printf("Running Parallel MTTKRP Benchmark for %s.\n", tensor_file);
    } else if (alg == -2) {
        selected_mttkrp_func = mttkrp_serial;
        printf("Running Serial MTTKRP Benchmark %s.\n", tensor_file);
    } else {
        fprintf(stderr, "Invalid algorithm value: %d. Expected -2 or -1.\n", alg);
        CU_cleanup_registry();
        return;
    }

	printf("Rank: %d\n", rank);
	printf("Threads: %d\n", num_threads);
	if (target_mode == -1 ) 
		printf("Target mode: all\n"); 
	else 
		printf("Mode: %d\n", target_mode);
	printf("Iterations: %d (skipping first warm-up)\n", num_iterations);
	printf("--------------------------------------------\n");

    if (target_mode != -1) {
        /* -------- Single mode benchmark -------- */
        double total_time = 0.0;

        for (int it = 0; it < num_iterations; ++it) {
            struct timespec start, end;
            clock_gettime(CLOCK_MONOTONIC, &start);

            matrix_t *computed = selected_mttkrp_func(global_tensor, global_factors, target_mode);

            clock_gettime(CLOCK_MONOTONIC, &end);
            double duration = (end.tv_sec - start.tv_sec) +
                              (end.tv_nsec - start.tv_nsec) / 1e9;

            printf("Mode %d Iteration %d Time: %.9f seconds\n", target_mode, it, duration);

            if (it > 0) total_time += duration;  // skip warmup
            free_matrix(computed);
        }

        double avg_time = total_time / (num_iterations - 1);
        printf("Mode %d MTTKRP Avg Time (excluding warm-up): %.9f seconds\n",
               target_mode, avg_time);

    } else {
        /* -------- All modes benchmark -------- */
        double grand_total = 0.0;

        for (int mode = 0; mode < global_tensor->ndims; ++mode) {
            double total_time_mode = 0.0;

            for (int it = 0; it < num_iterations; ++it) {
                struct timespec start, end;
                clock_gettime(CLOCK_MONOTONIC, &start);

                matrix_t *computed = selected_mttkrp_func(global_tensor, global_factors, mode);

                clock_gettime(CLOCK_MONOTONIC, &end);
                double duration = (end.tv_sec - start.tv_sec) +
                                  (end.tv_nsec - start.tv_nsec) / 1e9;

                printf("Mode %d Iteration %d Time: %.9f seconds\n", mode, it, duration);

                if (it > 0) total_time_mode += duration;  // skip warmup
                free_matrix(computed);
            }

            double avg_time_mode = total_time_mode / (num_iterations - 1);
            grand_total += avg_time_mode;
            printf("Mode %d MTTKRP Avg Time (excluding warm-up): %.9f seconds\n",
                   mode, avg_time_mode);
        }

        double overall_avg = grand_total / global_tensor->ndims;
        printf("Overall Average MTTKRP Time across %d modes: %.9f seconds\n",
               global_tensor->ndims, overall_avg);
    }

    suite_cleanup();
    CU_cleanup_registry();
}

/* Suite initialization: read all input files */
int suite_bench_init(const char *tensor_filename, int zero_base, int rank, int nnz) {

    // Read tensor
    FILE *file = fopen(tensor_filename, "r");
    if (!file) {
        perror("Error opening tensor file");
        exit(1);
    }
    global_tensor = read_tensor_file_with_base_fast(file, zero_base, nnz);
    fclose(file);
    if (!global_tensor) return 1;

    /* Allocate and generate factor matrices*/
    global_matrix_count = global_tensor->ndims;
    global_factors = (matrix_t **) MALLOC(sizeof(matrix_t *) * global_matrix_count);
    
    if (!global_factors) {
        fprintf(stderr, "Error allocating factor matrices\n");
        return 1;
    }

    for (int i = 0; i < global_matrix_count; ++i) {
        size_t rows = global_tensor->dims[i];
        global_factors[i] = new_random_matrix(rows, rank, 0.0, 1.0);  // example range [0,1]
        if (!global_factors[i]) {
            fprintf(stderr, "Error generating factor matrix for mode %d\n", i);
            return 1;
        }
    }

    return 0;
}


/* Verify if computed answers are equal to matlab's */
void CUnit_verify_mttkrp(const char *tensor_file, const char *factor_file, const char *mttkrp_file, int alg,int zero_base,int nnz) {
    CU_initialize_registry();

    if (suite_verify_init(tensor_file, factor_file, mttkrp_file,zero_base,nnz)) {
        fprintf(stderr, "Suite initialization failed.\n");
        CU_cleanup_registry();
        return;
    }

    CU_pSuite pSuite = CU_add_suite("MTTKRP Test", 0, 0);

    if (alg == -2) {
        selected_mttkrp_func = mttkrp_serial;
        printf("Running Serial MTTKRP Test\n");
    } else if (alg == -1) {
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

/* Compute MTTKRP over all modes */
matrix_t **get_mttkrp_results(struct hacoo_tensor *t, matrix_t **factor_matrices, int matrix_count, mttkrp_func_t f) {
    matrix_t **results = (matrix_t **)MALLOC(sizeof(matrix_t *) * t->ndims);
    for (int i = 0; i < matrix_count; i++) {
        results[i] = f(t, factor_matrices, i);
    }
    return results;
}

/* Suite initialization for verify MTTKRP: read all input files */
int suite_verify_init(const char *tensor_filename, const char *factor_filename, const char *mttkrp_filename, int zero_base,int nnz) {

    // Read tensor
    FILE *file = fopen(tensor_filename, "r");
    if (!file) {
        perror("Error opening tensor file");
        exit(1);
    }
    global_tensor = read_tensor_file_with_base_fast(file, zero_base, nnz);
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
