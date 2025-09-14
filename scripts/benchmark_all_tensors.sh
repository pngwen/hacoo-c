#!/bin/bash

#SBATCH --job-name=run-all-tensors
#SBATCH --partition=cpu
#SBATCH --output=all_tensors_output.txt
#SBATCH --error=all_tensors_error.txt
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128

# Tensors to run
tensors=("uber")
#tensors=("uber" "chicago" "lbnl" "nips" "nell-2" "enron")

# Path to benchmark script (must be executable and NOT include SBATCH headers)
benchmark_script="$HOME/haccoo-c/scripts/hacoo_mttkrp_bench.sh"

for tensor in "${tensors[@]}"; do
    echo "Running benchmark for tensor: $tensor"
    TENSOR_NAME="$tensor" bash "$benchmark_script"
done
