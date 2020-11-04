import numpy as np
import scqtl_tf
import scqtl_simulation
import pdb
# Generate some ZINB-distributed counts
num_samples = 1000
umi = np.concatenate([scqtl_simulation.simulate(
  num_samples=num_samples,
  size=1e5,
  seed=trial)[0][:,:1] for trial in range(10)], axis=1)
size_factor = 1e5 * np.ones((num_samples, 1))

# Generate a null design matrix
design = np.zeros((num_samples, 1))

# Map all samples to one individual/condition, i.e. one set of ZINB parameters
onehot = np.ones((num_samples, 1))

# Find the NB MLE
# Important: casting to float32 is required
init = scqtl_tf.fit(
  umi=umi.astype(np.float32),
  onehot=onehot.astype(np.float32),
  design=design.astype(np.float32),
  size_factor=size_factor.astype(np.float32),
  learning_rate=1e-3,
  max_epochs=20000,
  verbose=True,
)

# Find the ZINB MLE, starting from the NB MLE
log_mu, log_phi, logodds, nb_llik, zinb_llik = scqtl_tf.fit(
  umi=umi.astype(np.float32),
  onehot=onehot.astype(np.float32),
  design=design.astype(np.float32),
  size_factor=size_factor.astype(np.float32),
  learning_rate=1e-3,
  max_epochs=20000,
  warm_start=init[:3],
  verbose=True)

