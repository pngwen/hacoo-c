#!/bin/bash

#SBATCH --job-name=uber_mttkrp         # Job name
#SBATCH --partition=cpu
#SBATCH --output=output.txt             # Standard output file
#SBATCH --error=error.txt               # Standard error file
#SBATCH --time=01:00:00                 # Time limit
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks-per-node=1             # One task (sequential execution)
#SBATCH --cpus-per-task=64              # Match to OMP_NUM_THREADS

cd /u/jcharles1/haccoo-c

# Set OpenMP thread count
export OMP_NUM_THREADS=64

# Number of iterations
num_iterations=5

# Arrays to store run times
serial_times=()
parallel_times=()

echo "Starting Serial MTTKRP Tests..."
for ((i=1; i<=num_iterations; i++)); do
    echo "Serial Iteration $i..."

    start_time=$(date +%s.%N)
    ./run.sh uber 0
    end_time=$(date +%s.%N)

    duration=$(echo "$end_time - $start_time" | bc)
    serial_times+=("$duration")
    echo "Serial Iteration $i took ${duration} seconds."
done

# Calculate average serial time
serial_total=0
for t in "${serial_times[@]}"; do
    serial_total=$(echo "$serial_total + $t" | bc)
done
serial_avg=$(echo "$serial_total / ${#serial_times[@]}" | bc -l)

echo ""
echo "Average Serial MTTKRP Time: ${serial_avg} seconds"
echo ""

echo "Starting Parallel MTTKRP Tests..."
for ((i=1; i<=num_iterations; i++)); do
    echo "Parallel Iteration $i..."

    start_time=$(date +%s.%N)
    ./run.sh uber 1
    end_time=$(date +%s.%N)

    duration=$(echo "$end_time - $start_time" | bc)
    parallel_times+=("$duration")
    echo "Parallel Iteration $i took ${duration} seconds."
done

# Calculate average parallel time
parallel_total=0
for t in "${parallel_times[@]}"; do
    parallel_total=$(echo "$parallel_total + $t" | bc)
done
parallel_avg=$(echo "$parallel_total / ${#parallel_times[@]}" | bc -l)

echo ""
echo "Average Parallel MTTKRP Time: ${parallel_avg} seconds"
echo ""

echo "All iterations finished!"
