import numpy as np 
import os
import sys
import pdb

def prepare_expression(expression_file, output_expression):
	f = open(expression_file)
	t = open(output_expression, 'w')
	for line in f:
		line = line.rstrip()
		t.write(line + '\n')
	f.close()
	t.close()

def prepare_meta_data_one_hot_encoding(meta_data_file, covariate_name, output_file):
	f = open(meta_data_file)
	head_count = 0
	cell_names = []
	covariate_names = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if len(data) != 94:
			print('assumption error')
			pdb.set_trace()
		if head_count == 0:
			head_count = head_count + 1
			column_index_arr = np.where(np.asarray(data) == covariate_name)[0]
			if len(column_index_arr) != 1:
				print('assumption error: covariate name does not match only 1 column')
				pdb.set_trace()
			column_index = column_index_arr[0]
			continue
		cell_name = data[0]
		cell_names.append(cell_name)
		covariate_name = data[column_index]
		covariate_names.append(covariate_name)
	f.close()
	cell_names = np.asarray(cell_names)
	covariate_names = np.asarray(covariate_names)

	if len(cell_names) != len(np.unique(cell_names)):
		print('assumption error')

	ordered_covariates = np.unique(covariate_names)
	mapping = {}
	for i, cov in enumerate(ordered_covariates):
		mapping[cov] = i
	# output handle
	t = open(output_file,'w')
	# print header
	t.write('cell_id\t' + '\t'.join(ordered_covariates) + '\n')
	# Stream cells
	for index, cell_name in enumerate(cell_names):
		covariate_name = covariate_names[index]
		one_hot_vec = np.zeros(len(ordered_covariates))
		one_hot_vec[mapping[covariate_name]] = 1.0
		t.write(cell_name + '\t' + '\t'.join(one_hot_vec.astype(str)) + '\n')
	t.close()

def prepare_size_factors(meta_data_file, output_size_factor_file):
	f = open(meta_data_file)
	head_count = 0
	t = open(output_size_factor_file, 'w')
	# print header
	t.write('cell_id\tsize_factor\n')
	# stream cells
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			print(data[71])
			continue
		# Error checking
		if len(data) != 94:
			print('assumption eroror')
			pdb.set_trace()
		size_factor = data[71]
		cell_id = data[0]
		t.write(cell_id + '\t' + size_factor + '\n')
	f.close()
	t.close()

def prepare_design_matrix_from_expression_pcs(pc_file, num_pcs, meta_data_file, output_design_matrix_file):
	# First extract cell ids
	cell_ids = []
	head_count = 0
	f = open(meta_data_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		cell_ids.append(data[0])
	f.close()
	cell_ids = np.asarray(cell_ids)
	f = open(pc_file)
	t = open(output_design_matrix_file, 'w')
	# Print header
	col_names = []
	for col_num in range(num_pcs):
		col_names.append('cov' + str(col_num))
	col_names = np.asarray(col_names)
	t.write('cell_id\t' + '\t'.join(col_names) + '\n')
	cell_counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		t.write(cell_ids[cell_counter] + '\t' + '\t'.join(data[:num_pcs]) + '\n')
		cell_counter = cell_counter + 1
	f.close()
	t.close()
	if len(cell_ids) != cell_counter:
		print('assumption errror: mismatch in number of cells between covariate file and pc file')
		pdb.set_trace()

expression_file = sys.argv[1]
meta_data_file = sys.argv[2]
pc_file = sys.argv[3]
visualize_processed_expression_dir = sys.argv[4]


version_name = 'first_pass'
output_root = visualize_processed_expression_dir + version_name + '_vOutlier_input_data_'

#################################
# Prepare gene expression data
##################################
output_expression = output_root + 'expression.txt'
#prepare_expression(expression_file, output_expression

#################################
# Prepare one-hot encoding individaul matrix
##################################
output_one_hot_encoding_file = output_root + 'donor_one_hot_encoding.txt'
#prepare_meta_data_one_hot_encoding(meta_data_file, 'donor_long_id', output_one_hot_encoding_file)

#################################
# Prepare size factors
##################################
output_size_factor_file = output_root + 'size_factors.txt'
#prepare_size_factors(meta_data_file, output_size_factor_file)

#################################
# Prepare design matrix from expression pcs
##################################
num_pcs = 3
output_design_matrix_file = output_root + 'design_matrix_from_' + str(num_pcs) + '_pcs.txt'
#prepare_design_matrix_from_expression_pcs(pc_file, num_pcs, meta_data_file, output_design_matrix_file)
