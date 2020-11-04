#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --mem=10GB
#SBATCH --partition=lrgmem
#SBATCH --nodes=1


log_cpm_counts_file="$1"
meta_data_file="$2"
processed_expression_dir="$3"
visualize_processed_expression_dir="$4"

python process_expression.py $log_cpm_counts_file $meta_data_file $processed_expression_dir



module load R/3.5.1
Rscript visualize_processed_expression.R $processed_expression_dir $visualize_processed_expression_dir

