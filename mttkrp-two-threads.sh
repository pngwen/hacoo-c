#!/bin/bash

#SBATCH --job-name=chicago_mttkrp_scaling
#SBATCH --partition=cpu
#SBATCH --output=output.txt
#SBATCH --error=error.txt
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128

cd /u/jcharles1/haccoo-c

tensor="chicago"
tensor_file="MTTKRP_test/${tensor}/${tensor}_processed.txt"
factor_file="MTTKRP_test/${tensor}/factor_matrices.txt"
expected_file="MTTKRP_test/${tensor}/mttkrp_answers.txt"

echo "Tensor: $tensor"
echo ""
echo "MTTKRP: parallel (2 threads, 1 iteration)"
echo -e "threads\tmode 1\tmode 2\tmode 3\tmode 4\tsum"

export OMP_NUM_THREADS=2
output=$(./mttkrp_test "$tensor_file" "$factor_file" "$expected_file" 1 | grep 'Mode [0-9]' | awk '{print $5}')
readarray -t modes <<< "$output"

printf "2"
total=0
for ((j=0; j<${#modes[@]}; j++)); do
    printf "\t%.2f" "${modes[j]}"
    total=$(echo "$total + ${modes[j]}" | bc)
done
printf "\t%.2f\n" "$total"
