#!/bin/bash -l

#SBATCH
#SBATCH --time=12:00:00
#SBATCH --partition=gpuk80
#SBATCH --cpus-per-task=6
#SBATCH --gres=gpu:2
#SBATCH --ntasks-per-node=2




module load tensorflow/1.10.1-gpu-py3



input_expression_file="$1"
input_size_factor_file="$2"
input_one_hot_file="$3"
input_design_matrix_file="$4"
output_root="$5"


python vOutlier.py $input_expression_file $input_size_factor_file $input_one_hot_file $input_design_matrix_file $output_root
