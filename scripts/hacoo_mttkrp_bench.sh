#!/bin/bash

# Script to be run on a list of tensors using benchmark_all_tensors.slurm

tensor="${TENSOR_NAME:-test2}"  # default to test2 if not set
tensor_file="$HOME/haccoo-c/tensors/${tensor}.tns"
binary="$HOME/haccoo-c/hacoo_mttkrp"
rank=10
num_iterations=6

# Create timestamped log directory
timestamp=$(date +%Y%m%d_%H%M%S)
log_dir="$HOME/haccoo-c/benchmarking_logs/${tensor}_r${rank}_$timestamp"
mkdir -p "$log_dir"
results_file="$log_dir/results_${tensor}_r${rank}_$timestamp.tsv"

echo -e "Tensor: $tensor\n"
echo "Saving raw logs to $log_dir"
echo "Timing results will be written to $results_file"

# Determine number of modes
export OMP_NUM_THREADS=1
cmd="$binary -i \"$tensor_file\" -b -r $rank -t 1"
echo "Running: $cmd"
eval $cmd > "$log_dir/serial_probe.txt" 2>&1
output=$(awk '/Mode [0-9]+ Iteration/ { print $(NF-1) }' "$log_dir/serial_probe.txt")
readarray -t temp_modes <<< "$output"
num_modes=${#temp_modes[@]}

# Write header to TSV file
{
  echo "# Tensor: $tensor"
  echo "# Rank: $rank"
  header="threads\titeration"
  for ((j=0; j<num_modes; j++)); do
      header+="\tmode $j"
  done
  echo -e "$header"
} > "$results_file"

# PARALLEL TEST
echo "Running parallel MTTKRP..."

for threads in 8 16 32 64 128; do
    export OMP_NUM_THREADS=$threads
    log_file="$log_dir/parallel_t${threads}.txt"
    cmd="$binary -i \"$tensor_file\" -b -a -1 -r $rank -t $threads -n $num_iterations"
    echo "Running: $cmd"
    eval $cmd > "$log_file" 2>&1
    output=$(awk '/Mode [0-9]+ Iteration/ { print $(NF-1) }' "$log_file")
    readarray -t modes <<< "$output"

    iter=1
    for ((offset=0; offset<${#modes[@]}; offset+=num_modes)); do
        printf "%d\titeration %d" "$threads" "$iter" >> "$results_file"
        for ((j=0; j<num_modes; j++)); do
            val=$(printf "%.6f" "${modes[offset+j]}")
            printf "\t%.6f" "$val" >> "$results_file"
        done
        echo >> "$results_file"
        ((iter++))
    done
done
