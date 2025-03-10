#!/bin/bash

#SBATCH --job-name=uber_mttkrp      # Job name
#SBATCH --partition=cpu
#SBATCH --output=output.txt        # Standard output file
#SBATCH --error=error.txt          # Standard error file
#SBATCH --time=01:00:00            # Time limit
#SBATCH --nodes=1                  # Number of nodes
#SBATCH --ntasks-per-node=1         # One task (sequential execution)
#SBATCH --cpus-per-task=132       # match to OMP_NUM_THREADS

# Set OpenMP thread count
export OMP_NUM_THREADS=132

# Run the script multiple times using a loop
for i in {1..5}; do
    echo "Running iteration $i..."
    ./run.sh uber 0 # Replace with your actual script
    echo "Iteration $i completed."
done

echo "All iterations finished!"

# Run the script multiple times using a loop
for i in {1..5}; do
    echo "Running iteration $i..."
    ./run.sh uber 1 # Replace with your actual script
    echo "Iteration $i completed."
done

echo "All iterations finished!"
