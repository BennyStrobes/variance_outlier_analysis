import numpy as np 
import os
import sys
import pdb


def create_mapping_from_cell_id_to_size_factors(meta_data_file):
	f = open(meta_data_file)
	head_count = 0
	aa = []
	bb = []
	mapping = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			print(data[70])
			#pdb.set_trace()
			continue
		# Error checking
		if len(data) != 94:
			print('assumption eroror')
			pdb.set_trace()
		#pdb.set_trace()
		#size_factor = float(data[68])
		#size_factor = (np.power(10.0, float(data[17])) -1)/(np.power(10.0,6))
		size_factor = float(data[71])
		cell_id = data[0]
		if cell_id in mapping:
			print('assumption error')
			pdb.set_trace()
		mapping[cell_id] = size_factor
	f.close()
	return mapping

def transform_log_cpm_counts_to_raw_counts(log_cpm_counts_file, cell_id_to_size_factors, raw_counts_expression_file):
	f = open(log_cpm_counts_file)
	t = open(raw_counts_expression_file, 'w')
	# Initialze ordered array of size factors
	size_factors = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			for i, cell_id in enumerate(data):
				size_factor = cell_id_to_size_factors[cell_id]
				size_factors.append(size_factor)
			size_factors = np.asarray(size_factors)
			t.write('gene_id\t' + '\t'.join(data) + '\n')
			continue
		# Error checking
		if len(data) != 36045:
			print('assumption error')
			pdb.set_trace()
		# Extract relevent fields
		gene_id = data[0]
		log_cpm = np.asarray(data[1:]).astype(float)
		cpm = np.power(2.0, log_cpm) - 1.0
		# Something slightly wrong here..
		counts = np.round(cpm*size_factors/np.mean(size_factors)).astype(int)
		t.write(gene_id + '\t' + '\t'.join(counts.astype(str)) + '\n')
	f.close()
	t.close()

def filter_expression_to_day_0_cells(meta_data_file, raw_counts_expression_file, raw_counts_day_0_expression_file, new_header):
	day_0_cells = {}
	f = open(meta_data_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		day = data[6]
		if day != 'day0':
			continue
		day_0_cells[data[0]] = 0
	f.close()
	f = open(raw_counts_expression_file)
	t = open(raw_counts_day_0_expression_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			valid_indices = []
			if new_header == False:
				for i, cell_id in enumerate(data[1:]):
					if cell_id in day_0_cells:
						valid_indices.append(i)
				valid_indices = np.asarray(valid_indices)
				t.write(data[0] + '\t' + '\t'.join(np.asarray(data[1:])[valid_indices]) + '\n')
			else:
				for i, cell_id in enumerate(data):
					if cell_id in day_0_cells:
						valid_indices.append(i)
				valid_indices = np.asarray(valid_indices)
				t.write('gene_id' + '\t' + '\t'.join(np.asarray(data)[valid_indices]) + '\n')
			continue
		t.write(data[0] + '\t' + '\t'.join(np.asarray(data[1:])[valid_indices]) + '\n')
	f.close()
	t.close()

def filter_meta_data_file_to_match_expression_file(raw_counts_day_0_expression_file, meta_data_file, meta_data_day_0_file):
	cell_to_meta_data_mapping = {}
	f = open(meta_data_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			header = data
			continue
		cell_id = data[0]
		cell_to_meta_data_mapping[cell_id] = data
	f.close()
	ordered_cells = []
	f = open(raw_counts_day_0_expression_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		ordered_cells = data[1:]
		break
	f.close()
	t = open(meta_data_day_0_file,'w')
	t.write('cell_id\t' + '\t'.join(header) + '\n')
	for cell in ordered_cells:
		meta_data = cell_to_meta_data_mapping[cell]
		t.write('\t'.join(meta_data) + '\n')
	t.close()

def run_pca(expression_file, num_pcs, output_loadings_file, output_ve_file):
	# Load in data
	X_full = np.loadtxt(expression_file, dtype=str, delimiter='\t', comments='*')
	X = X_full[1:,1:].astype(float)
	# Standardize each gene
	num_genes = X.shape[0]
	for gene_num in range(num_genes):
		X[gene_num,:] = (X[gene_num,:] - np.mean(X[gene_num,:]))/np.std(X[gene_num,:])

	# Run PCA (via SVD)
	uuu, sss, vh = np.linalg.svd(X, full_matrices=False)
	svd_loadings = np.transpose(vh)[:,:num_pcs]

	# Save to output file
	np.savetxt(output_loadings_file, svd_loadings, fmt="%s", delimiter='\t')

	# Compute variance explained
	ve = (np.square(sss)/np.sum(np.square(sss)))[:num_pcs]
	np.savetxt(output_ve_file, ve, fmt="%s", delimiter='\n')

log_cpm_counts_file = sys.argv[1]
meta_data_file = sys.argv[2]
processed_expression_dir = sys.argv[3]


# First create mapping from cell id to size factors
cell_id_to_size_factors = create_mapping_from_cell_id_to_size_factors(meta_data_file)

# Transform log(cpm+1) expression matrix to raw counts expression matrix
raw_counts_expression_file = processed_expression_dir + 'raw_counts.txt'
transform_log_cpm_counts_to_raw_counts(log_cpm_counts_file, cell_id_to_size_factors, raw_counts_expression_file)

log_cpm_counts_day_0_expression_file = processed_expression_dir + 'log_cpm_counts_day_0.txt'
filter_expression_to_day_0_cells(meta_data_file, log_cpm_counts_file, log_cpm_counts_day_0_expression_file, True)

raw_counts_day_0_expression_file = processed_expression_dir + 'raw_counts_day_0.txt'
filter_expression_to_day_0_cells(meta_data_file, raw_counts_expression_file, raw_counts_day_0_expression_file, False)

meta_data_day_0_file = processed_expression_dir + 'cell_meta_data_day_0.txt'
filter_meta_data_file_to_match_expression_file(raw_counts_day_0_expression_file, meta_data_file, meta_data_day_0_file)

pca_loadings_day_0_file = processed_expression_dir + 'pca_loadings_day_0.txt'
pca_ve_day_0_file = processed_expression_dir + 'pca_variance_explained_day_0.txt'
num_pcs = 10
run_pca(log_cpm_counts_day_0_expression_file, num_pcs, pca_loadings_day_0_file, pca_ve_day_0_file)






