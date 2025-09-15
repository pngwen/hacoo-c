#!/bin/bash

#SBATCH --job-name=run-all-tensors
#SBATCH --partition=cpu
#SBATCH --output=all_tensors_output.txt
#SBATCH --error=all_tensors_error.txt
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128

# Default values
tensor="test3"
nnz="24"
#tensor="uber"
#nnz="3309490"

# Path to benchmark script (must be executable and NOT include SBATCH headers)
benchmark_script="$HOME/haccoo-c/scripts/new_hacoo_mttkrp_bench.sh"

echo "Running benchmark for tensor: $tensor"
if [ -n "$nnz" ]; then
    echo "Using number of nonzeros: $nnz"
    TENSOR_NAME="$tensor" NUM_NONZEROS="$nnz" bash "$benchmark_script"
else
    TENSOR_NAME="$tensor" bash "$benchmark_script"
fi
