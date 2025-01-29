#!/bin/bash

#echo Hello World

#   Run main with the following arguments:
#   global_argv[1] - name of file to read factor matrices from
#   global_argv[2] - name of file to read mttkrp results from

#Folder: mttkrp_test1
#./hacoo_test mttkrp_test1/factor_matrices.txt mttkrp_test1/mttkrp_answers.txt < mttkrp_test1/sptensor_data.tns

#Folder: mttkrp_test2
#./hacoo_test mttkrp_test2/factor_matrices.txt mttkrp_test2/mttkrp_answers.txt < mttkrp_test2/sptensor_data.tns

#Folder: uber_test
./hacoo_test uber_test/factor_matrices.txt uber_test/mttkrp_answers.txt < uber_test/uber.tns
