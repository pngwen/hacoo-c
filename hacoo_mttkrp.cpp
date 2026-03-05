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

/* Print usage guide */
void print_usage(const char *progname);

/* CUnit suite initialization and cleanup */
int suite_bench_init(const char *tensor_filename, const char *extension, int rank);
int suite_cleanup();
int generate_factor_matrices();

void CUnit_mttkrp_bench(const char *tensor_file, const char *extension, const char *output_file, int alg, int target_mode, 
                        int rank, int num_threads, int num_iterations, int nnz, int run_bench);

/* Globals */
struct hacoo_tensor *global_tensor = NULL;
matrix_t **global_factors = NULL;
int global_matrix_count = 0;

void print_usage(const char *progname) {
    printf("Usage: %s [OPTIONS]\n", progname);
    printf("Options:\n");
    printf("  -i or --input          Input tensor file (.tns or .hacoo); assumes indexes are 1-based\n");
    printf("  -o or --output         Output file name\n");
    printf("  -m or --mode           Target mode (-1:loop all modes, default; or specify a mode, e.g., 0 or 1 or 2 for third-order tensors.))\n");
    printf("  -s or --itrs           Number of iterations (1:default)\n");
    printf("  -d or --dev-id         MTTKRP algorithm (-2:sequential, default; -1:OpenMP parallel)\n");
    printf("  -r or --rank           Number of matrix columns (16:default)\n");
    printf("  -b or --bench          Run benchmark mode\n");
    printf("  -h or --help           Display this help message\n");
    printf("OpenMP options:\n");
    printf("  -t or --nt Number of threads (1:default)      \n");
    printf("\n");
}

/* Main function */
int main(int argc, char *argv[]) {
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

    char *tensor_file = NULL;
    char *output_file = NULL;
    char extension[10] = "tns"; // default COO file type
    int dev_id = -2; //default sequential
    int rank = 16;
    int target_mode = -1; //default all modes
    int run_bench = 0;
    int num_threads = 1;
	int num_iterations = 1;
	int nnz = 128;

    int opt;
	const char* const short_opt = "i:o:m:s:d:r:bht:";
	static struct option long_options[] = {
        {"input",       required_argument, 0, 'i'},
        {"output",      required_argument, 0, 'o'},
        {"mode",        required_argument, 0, 'm'},
        {"itrs",        required_argument, 0, 's'},
        {"dev-id",      required_argument, 0, 'd'},
		{"rank",        required_argument, 0, 'r'},
		{"bench",       no_argument,       0, 'b'},
        {"help",        no_argument,       0, 'h'},
        {"nt",          required_argument, 0, 't'},
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
            case 'o':
                output_file = optarg;
                break;
            case 'm':
                target_mode = atoi(optarg);
                break;    
			case 's':
				num_iterations = atoi(optarg);
				if (num_iterations <= 0) {
					fprintf(stderr, "Invalid number of iterations: %s\n", optarg);
					exit(1);
				}
				break;
            case 'r':
                rank = atoi(optarg);
                break;
            case 'd':
                dev_id = atoi(optarg);
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

    //determine file type (whatever's after the .)
    const char *ext = strrchr(tensor_file, '.');
    if (ext != NULL && ext != tensor_file) {
        strncpy(extension, ext + 1, sizeof(extension) - 1);
        extension[sizeof(extension) - 1] = '\0';
    }

    omp_set_num_threads(num_threads);
    openblas_set_num_threads(num_threads);
    CUnit_mttkrp_bench(tensor_file, extension, output_file, dev_id, target_mode, rank, num_threads, num_iterations, nnz, run_bench);

    return 0;
}

void CUnit_mttkrp_bench(const char *tensor_file, const char *extension,  const char *output_file, int alg,
                        int target_mode, int rank, int num_threads, int num_iterations, int nnz, int run_bench) {
    // Initialize CUnit
    CU_initialize_registry();
    
    if (suite_bench_init(tensor_file, extension, rank)) {
        fprintf(stderr, "Suite initialization failed.\n");
        CU_cleanup_registry();
        return;
    }

    // Select MTTKRP implementation
    if (alg == -1) {
        selected_mttkrp_func = hacoo_mttkrp;
        printf("Running Parallel MTTKRP for %s.\n", tensor_file);
    } else if (alg == -2) {
        selected_mttkrp_func = hacoo_mttkrp_serial;
        printf("Running Serial MTTKRP for %s.\n", tensor_file);
    } else {
        fprintf(stderr, "Invalid value: %d. Expected -2 or -1.\n", alg);
        CU_cleanup_registry();
        return;
    }

	printf("Rank: %d\n", rank);
	printf("Threads: %d\n", num_threads);
	if (target_mode == -1 ) { printf("Target mode: all\n"); }
    else { printf("Mode: %d\n", target_mode); }

    if(run_bench) { printf("Iterations: %d (skipping first warm-up)\n", num_iterations); }
	printf("--------------------------------------------\n");

    if(run_bench) { num_iterations++; }  // Increment to account for warm-up
    struct timespec start, end;

    if (target_mode != -1) {
        /* -------- Single mode -------- */
        double total_time = 0.0;

        for (int it = 0; it < num_iterations; ++it) {
            clock_gettime(CLOCK_MONOTONIC, &start);

            matrix_t *computed = selected_mttkrp_func(global_tensor, global_factors, target_mode);

            if(run_bench) {
                clock_gettime(CLOCK_MONOTONIC, &end);
                double duration = (end.tv_sec - start.tv_sec) +
                                (end.tv_nsec - start.tv_nsec) / 1e9;

                printf("Mode %d Iteration %d Time: %.9f seconds\n", target_mode, it, duration);
                if (it > 0) total_time += duration;  // skip warmup
            } else {
                // write output to file
                if (output_file) {
                    char filename[256];
                    snprintf(filename, sizeof(filename), "%s_mode_%d.txt", output_file, target_mode);
                    write_matrix_to_file(filename, computed);
                }
            }
            
            free_matrix(computed);
        }

        if(run_bench) {
            double avg_time = total_time / (num_iterations);
            printf("Mode %d MTTKRP Avg Time (excluding warm-up): %.9f seconds\n",
               target_mode, avg_time);
        }
    } else {
        /* -------- All modes -------- */
        double grand_total = 0.0;

        for (int mode = 0; mode < global_tensor->ndims; ++mode) {
            double total_time_mode = 0.0;

            for (int it = 0; it < num_iterations; ++it) {

                clock_gettime(CLOCK_MONOTONIC, &start);

                matrix_t *computed = selected_mttkrp_func(global_tensor, global_factors, mode);

                if(run_bench) {
                    clock_gettime(CLOCK_MONOTONIC, &end);
                    double duration = (end.tv_sec - start.tv_sec) +
                                  (end.tv_nsec - start.tv_nsec) / 1e9;

                    printf("Mode %d Iteration %d Time: %.9f seconds\n", mode, it, duration);
                    if (it > 0) total_time_mode += duration;  // skip warmup
                } else {
                    // write output to file
                    if (output_file) {
                        char filename[256];
                        snprintf(filename, sizeof(filename), "%s_mode_%d.txt", output_file, mode);
                        write_matrix_to_file(filename, computed);
                    }
                }
                
                free_matrix(computed);
            }

            if(run_bench) {
                double avg_time_mode = total_time_mode / (num_iterations-1);
                grand_total += avg_time_mode;
                printf("Mode %d MTTKRP Avg Time (excluding warm-up): %.9f seconds\n",
                    mode, avg_time_mode);
            }
        }

        if(run_bench) {
            double overall_avg = grand_total / global_tensor->ndims;
            printf("Overall Average MTTKRP Time across %lu modes: %.9f seconds\n",
               global_tensor->ndims, overall_avg);
        }
    }

    suite_cleanup();
    CU_cleanup_registry();
}

/* Suite initialization: read all input files */
int suite_bench_init(const char *tensor_filename, const char *extension, int rank) {

    // Read tensor
    FILE *file = fopen(tensor_filename, "r");
    if (!file) {
        perror("Error opening tensor file");
        exit(1);
    }

    if (strcmp(extension, "tns") == 0) {
        global_tensor = hacoo_read_tensor_file(file);
    } else { //assume extension .hacoo
        global_tensor = hacoo_read_htensor_file(file);
    }
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

/* Suite cleanup: free all loaded data */
int suite_cleanup() {
    if (global_tensor) {
        hacoo_free(global_tensor);
        global_tensor = NULL;
    }
    if (global_factors) {
        free_matrices(global_factors, global_matrix_count);
        global_factors = NULL;
    }
    global_matrix_count = 0;

    return 0;
}
