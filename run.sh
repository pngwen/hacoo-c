#!/bin/bash

# usage: ./run.sh uber 0 to run hacoo_test using serial mttkrp algorithm
# args:
# 1) name of tensor
# 2) mttkrp mode: 0 for serial, 1 parallel

# Check if the required arguments are provided
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <folder_name> <mttkrp_mode>"
    exit 1
fi

# Define folder and tensor names
folder_name="$1"
tensor_name="${folder_name}_processed.txt"

# Define variables for input files and execution mode
data_file="MTTKRP_test/$folder_name/$tensor_name"
factor_matrices_file="MTTKRP_test/$folder_name/factor_matrices.txt"
mttkrp_results_file="MTTKRP_test/$folder_name/mttkrp_answers.txt"
mttkrp_mode="$2" # 0 for serial, 1 for parallel

# Print what will run
echo "Running hacoo_test with:"
echo "Tensor file:         $data_file"
echo "Factor matrices:     $factor_matrices_file"
echo "Expected MTTKRP:     $mttkrp_results_file"
echo "MTTKRP mode:         $mttkrp_mode (0=serial, 1=parallel)"
echo ""

# Execute the hacoo_test program
./hacoo_test "$data_file" "$factor_matrices_file" "$mttkrp_results_file" "$mttkrp_mode"