#!/bin/bash

#SBATCH --job-name=test1_mttkrp_scaling
#SBATCH --partition=cpu
#SBATCH --output=output.txt
#SBATCH --error=error.txt
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128

cd /u/jcharles1/haccoo-c

tensor="test1"
tensor_file="MTTKRP_test/${tensor}/${tensor}_processed.txt"
factor_file="MTTKRP_test/${tensor}/factor_matrices.txt"
expected_file="MTTKRP_test/${tensor}/mttkrp_answers.txt"

num_iterations=5

echo "Tensor: $tensor"
echo ""

# SERIAL TESTS
echo "MTTKRP: serial"
echo -e "iteration\tmode 1\tmode 2\tmode 3\tmode 4\tsum"

# Initialize serial totals
serial_totals=(0 0 0 0)

for ((i=1; i<=num_iterations; i++)); do
    export OMP_NUM_THREADS=1
    output=$(./mttkrp_test "$tensor_file" "$factor_file" "$expected_file" 0 | grep 'Mode [0-9]' | awk '{print $5}')
    readarray -t modes <<< "$output"

    printf "iteration %d" "$i"
    total=0
    for ((j=0; j<${#modes[@]}; j++)); do
        printf "\t%.2f" "${modes[j]}"
        serial_totals[$j]=$(echo "${serial_totals[$j]} + ${modes[j]}" | bc)
        total=$(echo "$total + ${modes[j]}" | bc)
    done
    printf "\t%.2f\n" "$total"
done

# Print serial averages
echo -n "average"
total=0
for ((j=0; j<${#serial_totals[@]}; j++)); do
    avg=$(echo "${serial_totals[$j]} / $num_iterations" | bc -l)
    printf "\t%.2f" "$avg"
    total=$(echo "$total + $avg" | bc)
done
printf "\t%.2f\n" "$total"

echo ""
echo "MTTKRP: parallel"
echo -e "threads\titeration\tmode 1\tmode 2\tmode 3\tmode 4\tsum"

for threads in 1 2 4 8 16 32 64 128; do
    export OMP_NUM_THREADS="$threads"
    # Reset parallel totals
    parallel_totals=(0 0 0 0)

    for ((i=1; i<=num_iterations; i++)); do
        output=$(./mttkrp_test "$tensor_file" "$factor_file" "$expected_file" 1 | grep 'Mode [0-9]' | awk '{print $5}')
        readarray -t modes <<< "$output"

        printf "%d\titeration %d" "$threads" "$i"
        total=0
        for ((j=0; j<${#modes[@]}; j++)); do
            printf "\t%.2f" "${modes[j]}"
            parallel_totals[$j]=$(echo "${parallel_totals[$j]} + ${modes[j]}" | bc)
            total=$(echo "$total + ${modes[j]}" | bc)
        done
        printf "\t%.2f\n" "$total"
    done

    # Print averages for this thread count
    printf "%d\taverage" "$threads"
    total=0
    for ((j=0; j<${#parallel_totals[@]}; j++)); do
        avg=$(echo "${parallel_totals[$j]} / $num_iterations" | bc -l)
        printf "\t%.2f" "$avg"
        total=$(echo "$total + $avg" | bc)
    done
    printf "\t%.2f\n" "$total"
done
