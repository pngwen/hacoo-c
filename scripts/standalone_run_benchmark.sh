#!/bin/bash

tensor="uber"
tensor_file="$HOME/haccoo-c/tensors/${tensor}.tns"
binary="$HOME/haccoo-c/hacoo_mttkrp"
rank=10
num_iterations=3

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
output=$(awk '/Mode [0-9]+ MTTKRP Time:/ { print $(NF-1) }' "$log_dir/serial_probe.txt")
readarray -t temp_modes <<< "$output"
num_modes=${#temp_modes[@]}

# Write header to TSV file
{
  echo "# Tensor: $tensor"
  echo "# Rank: $rank"
  header="threads\titeration"
  for ((j=1; j<=num_modes; j++)); do
      header+="\tmode $j"
  done
  echo -e "$header"
} > "$results_file"

# SERIAL TEST
echo "Running serial MTTKRP..."
for ((j=0; j<num_modes; j++)); do
    serial_totals[$j]=0
done

for ((i=1; i<=num_iterations; i++)); do
    export OMP_NUM_THREADS=1
    log_file="$log_dir/serial_iter${i}.txt"
    cmd="$binary -i \"$tensor_file\" -b -r $rank -t 1"
    echo "Running: $cmd"
    eval $cmd > "$log_file" 2>&1
    output=$(awk '/Mode [0-9]+ MTTKRP Time:/ { print $(NF-1) }' "$log_file")
    readarray -t modes <<< "$output"

    printf "1\titeration %d" "$i" >> "$results_file"
    for ((j=0; j<num_modes; j++)); do
        val=$(printf "%.2f" "${modes[j]}")
        printf "\t%.2f" "$val" >> "$results_file"
        serial_totals[$j]=$(echo "${serial_totals[$j]} + $val" | bc)
    done
    echo >> "$results_file"
done

# Serial average
printf "1\taverage" >> "$results_file"
for ((j=0; j<num_modes; j++)); do
    avg=$(echo "scale=3; ${serial_totals[$j]} / $num_iterations" | bc)
    printf "\t%.3f" "$avg" >> "$results_file"
done
echo >> "$results_file"

# PARALLEL TEST
echo "Running parallel MTTKRP..."

for threads in 1; do
#for threads in 1 2 4 8 16 32 64 128; do
    export OMP_NUM_THREADS=$threads

    for ((j=0; j<num_modes; j++)); do
        parallel_totals[$j]=0
    done

    for ((i=1; i<=num_iterations; i++)); do
        log_file="$log_dir/parallel_t${threads}_iter${i}.txt"
        cmd="$binary -i \"$tensor_file\" -b -a -1 -r $rank -t $threads"
        echo "Running: $cmd"
        eval $cmd > "$log_file" 2>&1
        output=$(awk '/Mode [0-9]+ MTTKRP Time:/ { print $(NF-1) }' "$log_file")
        readarray -t modes <<< "$output"

        printf "%d\titeration %d" "$threads" "$i" >> "$results_file"
        for ((j=0; j<num_modes; j++)); do
            val=$(printf "%.2f" "${modes[j]}")
            printf "\t%.2f" "$val" >> "$results_file"
            parallel_totals[$j]=$(echo "${parallel_totals[$j]} + $val" | bc)
        done
        echo >> "$results_file"
    done

    printf "%d\taverage" "$threads" >> "$results_file"
    for ((j=0; j<num_modes; j++)); do
        avg=$(echo "scale=3; ${parallel_totals[$j]} / $num_iterations" | bc)
        printf "\t%.3f" "$avg" >> "$results_file"
    done
    echo >> "$results_file"
done
