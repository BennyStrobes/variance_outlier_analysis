##################
# Input Data
##################
log_cpm_counts_file="/work-zfs/abattle4/lab_data/sc_endo_diff/counts.tsv"
meta_data_file="/work-zfs/abattle4/lab_data/sc_endo_diff/cell_metadata_cols.tsv"



##################
# Output directories
##################
# Output root directory
output_root="/work-zfs/abattle4/bstrober/rare_variant_variance_outliers/preprocess_cuomo_data/"
# Directory containing processed expression data
processed_expression_dir=$output_root"processed_expression/"
# Directory containing visualizations of processed expression data
visualize_processed_expression_dir=$output_root"visualize_processed_expression/"
# Variance outlier input data dir
variance_outlier_input_data_dir=$output_root"variance_outlier_input_data/"
# Variance outlier input data dir
variance_outlier_results_dir=$output_root"variance_outlier_results/"


if false; then
sh preprocess_expression.sh $log_cpm_counts_file $meta_data_file $processed_expression_dir $visualize_processed_expression_dir
fi

expression_file=$processed_expression_dir"raw_counts_day_0.txt"
covariate_file=$processed_expression_dir"cell_meta_data_day_0.txt"
pc_file=$processed_expression_dir"pca_loadings_day_0.txt"
if false; then
sh prepare_variance_outlier_input_data.sh $expression_file $covariate_file $pc_file $variance_outlier_input_data_dir
fi





vOutlier_input_expression_file=$variance_outlier_input_data_dir"first_pass_vOutlier_input_data_expression.txt"
vOutlier_input_size_factor_file=$variance_outlier_input_data_dir"first_pass_vOutlier_input_data_size_factors.txt"
vOutlier_input_one_hot_file=$variance_outlier_input_data_dir"first_pass_vOutlier_input_data_donor_one_hot_encoding.txt"
vOutlier_input_design_matrix_file=$variance_outlier_input_data_dir"first_pass_vOutlier_input_data_design_matrix_from_3_pcs.txt"
vOutlier_output_root=$variance_outlier_results_dir"first_pass_vOutlier_results_"
sbatch run_vOutlier_analysis.sh $vOutlier_input_expression_file $vOutlier_input_size_factor_file $vOutlier_input_one_hot_file $vOutlier_input_design_matrix_file $vOutlier_output_root