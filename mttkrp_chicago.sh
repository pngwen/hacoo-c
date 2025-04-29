#!/bin/bash

#SBATCH --job-name=chicago_mttkrp         # Job name
#SBATCH --partition=cpu
#SBATCH --output=output.txt            # Standard output file
#SBATCH --error=error.txt              # Standard error file
#SBATCH --time=01:00:00                # Time limit
#SBATCH --nodes=1                      # Number of nodes
#SBATCH --ntasks-per-node=1            # One task (sequential execution)
#SBATCH --cpus-per-task=64             # Match to OMP_NUM_THREADS

cd /u/jcharles1/haccoo-c

# Set OpenMP thread count
export OMP_NUM_THREADS=64

# Number of iterations
num_iterations=5

echo "Starting Serial MTTKRP Tests..."
for ((i=1; i<=num_iterations; i++)); do
    echo "Serial Iteration $i..."
    ./run.sh chicago 0
done

echo ""
echo "Starting Parallel MTTKRP Tests..."
for ((i=1; i<=num_iterations; i++)); do
    echo "Parallel Iteration $i..."
    ./run.sh chicago 1
done

echo ""
echo "All iterations finished!"
