#!/bin/bash

#SBATCH --job-name=uber_mttkrp         # Job name
#SBATCH --partition=cpu
#SBATCH --output=output.txt            # Standard output file
#SBATCH --error=error.txt             # Standard error file
#SBATCH --time=01:00:00                # Time limit
#SBATCH --nodes=1                      # Number of nodes
#SBATCH --ntasks-per-node=1            # One task (sequential execution)
#SBATCH --cpus-per-task=128             # Match to OMP_NUM_THREADS

# Move to the directory containing the executable
cd /u/jcharles1/haccoo-c

# Number of iterations
num_iterations=5

# Function to run mttkrp_test
run_mttkrp_test () {
    folder_name="$1"
    mttkrp_alg="$2"
    tensor_name="${folder_name}_processed.txt"

    data_file="MTTKRP_test/$folder_name/$tensor_name"
    factor_matrices_file="MTTKRP_test/$folder_name/factor_matrices.txt"
    mttkrp_results_file="MTTKRP_test/$folder_name/mttkrp_answers.txt"

    echo "Running mttkrp_test with:"
    echo "Tensor file:         $data_file"
    echo "Factor matrices:     $factor_matrices_file"
    echo "Expected MTTKRP:     $mttkrp_results_file"
    echo "MTTKRP alg:          $mttkrp_alg (0=serial, 1=parallel)"
    echo ""

    ./mttkrp_test "$data_file" "$factor_matrices_file" "$mttkrp_results_file" "$mttkrp_alg"
}

echo "Starting Serial MTTKRP Tests..."
for ((i=1; i<=num_iterations; i++)); do
    echo "Serial Iteration $i..."
    run_mttkrp_test uber 0
done

echo ""
echo "Starting Parallel MTTKRP Tests..."

# Loop through powers of 2 for OMP_NUM_THREADS, up to 128
for threads in 1 2 4 8 16 32 64 128; do
    export OMP_NUM_THREADS=$threads
    echo "Running with $OMP_NUM_THREADS threads..."
    for ((i=1; i<=num_iterations; i++)); do
        echo "Parallel Iteration $i with $OMP_NUM_THREADS threads..."
        run_mttkrp_test uber 1
    done
done

echo ""
echo "All iterations finished!"
