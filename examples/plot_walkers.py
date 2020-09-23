import numpy as np
import h5py
import matplotlib.pyplot as plt

# path to mcmc chain file
ssp_dir = './'
mcmc = 'mcmc_chain.hdf5'

chain_file = h5py.File(ssp_dir+mcmc, 'r')

print('Parameters of walkers: ',chain_file.keys())

pars = ['logmass', 'logZ', 'tau', 'age', 'tau_V', 'logp', 'mag_u', 'mag_g', 'mag_r', 'mag_i', 'mag_z']
chain = {}
for par in pars:
    chain[par] = chain_file[par]

magnitudes = ['u', 'g', 'r', 'i', 'z']
mags = {'u': 22.30, 'g': 20.26, 'r': 18.79, 'i': 18.32, 'z': 18.04} #synthetic

nr_walkers = 50 #number of walkers

# Plot the walkers of stellar parameters + logp.
fig, axes = plt.subplots(ncols=1, nrows=len(chain)-len(mags))
fig.set_size_inches(12,9)

for j in np.arange(0, len(chain)-len(mags)):
	for i in range(nr_walkers):
		axes[j].plot(chain[pars[j]][i,:], color='black', alpha=0.3)
		axes[j].set_ylabel(pars[j])
		axes[j].set_xlim([0,300])

plt.savefig('./walkers_stel_parameters.pdf', overwrite=True)
plt.show()

# plot the walkers from the magnitudes.
fig2, axes2 = plt.subplots(ncols=1, nrows=len(mags))
fig2.set_size_inches(12,8)

for j in np.arange(len(chain)-len(mags), len(chain)):
	ax = j- len(mags)-2
	print(ax)
	for i in range(nr_walkers):
		axes2[ax].plot(chain[pars[j]][i,:], color='black', alpha=0.3)
		axes2[ax].set_ylabel(pars[j])
		axes2[ax].axhline(mags[magnitudes[ax]], color='r')
		axes2[ax].set_xlim([0,300])	

plt.savefig('./walkers_magnitudes.pdf', overwrite=True)
plt.show()

chain_file.close()

