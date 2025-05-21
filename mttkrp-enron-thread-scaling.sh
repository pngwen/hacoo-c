#!/bin/bash

#SBATCH --job-name=mttkrp-enron-thread-scaling
#SBATCH --partition=cpu
#SBATCH --output=output.txt
#SBATCH --error=error.txt
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128

cd /u/jcharles1/haccoo-c

tensor="enron"
tensor_file="MTTKRP_test/${tensor}/${tensor}_processed.txt"
factor_file="MTTKRP_test/${tensor}/factor_matrices.txt"
expected_file="MTTKRP_test/${tensor}/mttkrp_answers.txt"

num_iterations=5

echo -e "Tensor: $tensor\n"

# Determine number of modes from a single test run
export OMP_NUM_THREADS=1
first_output=$(./mttkrp_test "$tensor_file" "$factor_file" "$expected_file" 0 | grep 'Mode [0-9]' | awk '{print $5}')
readarray -t temp_modes <<< "$first_output"
num_modes=${#temp_modes[@]}

# Prepare header line
header="iteration"
for ((j=1; j<=num_modes; j++)); do
    header+="\tmode $j"
done

# SERIAL TEST
echo -e "MTTKRP: serial"
echo -e "$header"

# Initialize serial totals
for ((j=0; j<num_modes; j++)); do
    serial_totals[$j]=0
done

for ((i=1; i<=num_iterations; i++)); do
    export OMP_NUM_THREADS=1
    output=$(./mttkrp_test "$tensor_file" "$factor_file" "$expected_file" 0 | grep 'Mode [0-9]' | awk '{print $5}')
    readarray -t modes <<< "$output"

    printf "iteration %d" "$i"
    for ((j=0; j<num_modes; j++)); do
        val=$(printf "%.3f" "${modes[j]}")
        printf "\t%.2f" "$val"
        serial_totals[$j]=$(echo "${serial_totals[$j]} + ${val}" | bc)
    done
    echo
done

# Print serial averages
printf "average"
for ((j=0; j<num_modes; j++)); do
    avg=$(echo "scale=3; ${serial_totals[$j]} / $num_iterations" | bc)
    printf "\t%.3f" "$avg"
done
echo -e "\n"

# PARALLEL TEST
echo -e "MTTKRP: parallel"
echo -e "threads\t$header"

for threads in 1 2 4 8 16 32 64 128; do
    export OMP_NUM_THREADS="$threads"

    # Reset totals
    for ((j=0; j<num_modes; j++)); do
        parallel_totals[$j]=0
    done

    for ((i=1; i<=num_iterations; i++)); do
        output=$(./mttkrp_test "$tensor_file" "$factor_file" "$expected_file" 1 | grep 'Mode [0-9]' | awk '{print $5}')
        readarray -t modes <<< "$output"

        printf "%d\titeration %d" "$threads" "$i"
        for ((j=0; j<num_modes; j++)); do
            val=$(printf "%.3f" "${modes[j]}")
            printf "\t%.2f" "$val"
            parallel_totals[$j]=$(echo "${parallel_totals[$j]} + ${val}" | bc)
        done
        echo
    done

    printf "%d\taverage" "$threads"
    for ((j=0; j<num_modes; j++)); do
        avg=$(echo "scale=3; ${parallel_totals[$j]} / $num_iterations" | bc)
        printf "\t%.3f" "$avg"
    done
    echo
done
