#!/bin/bash
#SBATCH --job-name=loop_job        # Job name
#SBATCH --output=output.txt        # Standard output file
#SBATCH --error=error.txt          # Standard error file
#SBATCH --time=01:00:00            # Time limit
#SBATCH --partition=compute        # Partition name
#SBATCH --nodes=1                  # Number of nodes
#SBATCH --ntasks=1                 # One task (sequential execution)
#SBATCH --cpus-per-task=8          # OpenMP threads

# Load necessary modules
#module load my_software

# Set OpenMP thread count
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run the script multiple times using a loop
for i in {1..5}; do
    echo "Running iteration $i..."
    ./run.sh uber 0   # Replace with your actual script
    echo "Iteration $i completed."
done

echo "All iterations finished!"

# Run the script multiple times using a loop
for i in {1..5}; do
    echo "Running iteration $i..."
    ./run.sh uber 1   # Replace with your actual script
    echo "Iteration $i completed."
done

echo "All iterations finished!"
