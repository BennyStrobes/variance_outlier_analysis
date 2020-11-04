#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --mem=10GB
#SBATCH --partition=lrgmem
#SBATCH --nodes=1


expression_file="$1"
meta_data_file="$2"
pc_file="$3"
visualize_processed_expression_dir="$4"



python prepare_variance_outlier_input_data.py $expression_file $meta_data_file $pc_file $visualize_processed_expression_dir