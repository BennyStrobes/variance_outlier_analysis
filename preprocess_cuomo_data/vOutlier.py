import numpy as np
import scqtl_tf
import pdb
import sys




def load_input_data(input_size_factor_file, input_one_hot_file, input_design_matrix_file):
	# Size factors
	raw_size_factors = np.loadtxt(input_size_factor_file,dtype=str,delimiter='\t', comments='*')
	size_factors = np.reshape(raw_size_factors[1:,1], (-1,1)).astype(np.float32)
	size_factor_cell_names = raw_size_factors[1:,0]
	# One hot
	raw_one_hot = np.loadtxt(input_one_hot_file,dtype=str,delimiter='\t', comments='*')
	onehot = raw_one_hot[1:,1:].astype(np.float32)
	onehot_cell_names = raw_one_hot[1:,0]
	individual_id_names = raw_one_hot[0,1:]
	# Design matrix
	raw_design = np.loadtxt(input_design_matrix_file, dtype=str, delimiter='\t', comments='*')
	design = raw_design[1:,1:].astype(np.float32)
	design_cell_names = raw_design[1:,0]

	if np.array_equal(design_cell_names, onehot_cell_names) == False:
		print('assumption error: cell ids for design matrix and one hot do not correspond')
		pdb.set_trace()
	if np.array_equal(design_cell_names, size_factor_cell_names) == False:
		print('assumption error: cell ids for design matrix and design matrix do not correspond')
		pdb.set_trace()
	return size_factors, design, onehot, size_factor_cell_names, individual_id_names


######################################
# Run vOutlier analysis for single gene
#######################################
def vOutlier_for_single_gene(counts, size_factor, design, onehot):
	# Find the NB MLE
	init = scqtl_tf.fit(umi=counts, onehot=onehot, design=design, size_factor=size_factor, learning_rate=1e-3, max_epochs=20000, verbose=True)

	# Find the ZINB MLE, starting from the NB MLE
	log_mu, log_phi, logodds, nb_llik, zinb_llik = scqtl_tf.fit(umi=counts, onehot=onehot, design=design, size_factor=size_factor, learning_rate=1e-3, max_epochs=20000, warm_start=init[:3], verbose=True)

	mu = np.exp(log_mu)[:,0]
	phi = np.exp(log_phi)[:,0]
	logodds = logodds[:,0]
	return mu, phi, logodds

#####################
# Command line args
#####################
input_expression_file = sys.argv[1]
input_size_factor_file = sys.argv[2]
input_one_hot_file = sys.argv[3]
input_design_matrix_file = sys.argv[4]
output_root = sys.argv[5]


##############
# Load in input data
size_factor, design, onehot, cell_ids, individual_ids = load_input_data(input_size_factor_file, input_one_hot_file, input_design_matrix_file)

##############
# Open output file handles
t_mu = open(output_root + 'param_mu.txt', 'w')
t_phi = open(output_root + 'param_phi.txt', 'w')
t_logodds = open(output_root + 'param_logodds.txt','w')

t_mu.write('gene_id\t' + '\t'.join(individual_ids) + '\n')
t_phi.write('gene_id\t' + '\t'.join(individual_ids) + '\n')
t_logodds.write('gene_id\t' + '\t'.join(individual_ids) + '\n')


###############
# Loop through genes to run vOutlier analysis for each gene
f = open(input_expression_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		# HEADER
		head_count = head_count + 1
		if np.array_equal(data[1:], cell_ids) == False:
			print('assumption error: expression cell ids do not match covariate cell ids')
			pdb.set_trace()
		continue
	# Extract relevent data for thi sgene
	gene_id = data[0]
	print('##############################')
	print('vOutlier analysis for ' + gene_id)
	print('##############################')
	counts = np.reshape(np.asarray(data[1:]).astype(np.float32), (-1, 1))
	# Run vOutlier analysis for this gene
	mu, phi, logodds = vOutlier_for_single_gene(counts, size_factor, design, onehot)
	# print to output file handles
	t_mu.write(gene_id + '\t' + '\t'.join(mu.astype(str)) + '\n')
	t_phi.write(gene_id + '\t' + '\t'.join(phi.astype(str)) + '\n')
	t_logodds.write(gene_id + '\t' + '\t'.join(logodds.astype(str)) + '\n')
	# Flush
	t_mu.flush()
	t_phi.flush()
	t_logodds.flush()
f.close()
t_mu.close()
t_phi.close()
t_logodds.close()

