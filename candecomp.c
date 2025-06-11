#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cpd.h"
#include "hacoo.h"
#include "matrix.h"

#define DEFAULT_RANK 10
#define DEFAULT_MAX_ITER 1000

void print_usage(const char *program_name)
{
    printf("Usage: %s <filename> [--rank <rank>] [--max_iter <max_iter>]\n", program_name);
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        print_usage(argv[0]);
        return 1;
    }

    const char *filename = argv[1];
    unsigned int rank = DEFAULT_RANK;
    unsigned int max_iter = DEFAULT_MAX_ITER;

    // Parse optional arguments
    for (int i = 2; i < argc; i++)
    {
        if (strcmp(argv[i], "--rank") == 0 && i + 1 < argc)
        {
            rank = atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "--max_iter") == 0 && i + 1 < argc)
        {
            max_iter = atoi(argv[++i]);
        }
        else
        {
            print_usage(argv[0]);
            return 1;
        }
    }

    // Read the tensor from the file
    FILE *file = fopen(filename, "r");
    if(!file) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        return 1;
    }
    struct hacoo_tensor *tensor = read_tensor_file(file);
    if (tensor == NULL)
    {
        fprintf(stderr, "Error reading tensor from file: %s\n", filename);
        return 1;
    }

    // Perform CPD
    double tol = 1e-5;
    cpd_result_t *result= cpd(tensor, rank, max_iter, tol);

    // Print the factor matrices
    for (unsigned int i = 0; i < tensor->ndims; i++)
    {
        printf("Factor matrix %u:\n", i);
        print_matrix(result->factors[i]);
    }

    cpd_result_free(result);
    hacoo_free(tensor);

    return 0;
}
