#!/bin/bash

#echo Hello World

#   Run main with the following arguments:
#   global_argv[1] - name of file to read factor matrices from
#   global_argv[2] - name of file to read mttkrp results from

#Folder: mttkrp_test1
#./hacoo_test mttkrp_test1/factor_matrices.txt mttkrp_test1/mttkrp_answers.txt < mttkrp_test1/sptensor_data.tns

#Folder: mttkrp_test2
#./hacoo_test mttkrp_test2/factor_matrices.txt mttkrp_test2/mttkrp_answers.txt < mttkrp_test2/sptensor_data.tns

#Folder: mttkrp_test3 (1 based to 0-based calc MATLAB script, subtracts 1 from all indexes in hacoo.c library)
#./hacoo_test mttkrp_test3/factor_matrices.txt mttkrp_test3/mttkrp_answers.txt < mttkrp_test3/sptensor_data.tns

#Folder: mttkrp_test4 (rehash test)
#./hacoo_test mttkrp_test4/factor_matrices.txt mttkrp_test4/mttkrp_answers.txt < mttkrp_test4/sptensor_data.tns

#MAKE SURE YOU UNCOMMENT CODE THAT SUBTRACTS 1 FROM INDEXES
#DID YOU READ THAT
#
#
#Folder: mttkrp_test5 (index fixing test)
#./hacoo_test mttkrp_test5/factor_matrices.txt mttkrp_test5/mttkrp_answers.txt < mttkrp_test5/sptensor_data.tns

./hacoo_test uber_test/factor_matrices.txt uber_test/mttkrp_answers.txt < uber_test/uber.tns
