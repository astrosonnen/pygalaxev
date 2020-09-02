import numpy as np
import emcee
import ndinterp
import pygalaxev_cosmology
from scipy.interpolate import splrep
from scipy.stats import truncnorm
import h5py


griddir = '/data2/sonnenfeld/stellarpop_csp_models/bc03_basel_chabrier/'
work_dir = griddir

redshift = 0.1945

outname = 'mcmc_chain.hdf5'
# MCMC parameters
nstep = 500
nwalkers = 50
burnin = 200
thin = 1

bands = ['u', 'g', 'r', 'i']

# observed broad band magnitudes and uncertainties
mags = {'u': 20.45, 'g': 18.52, 'r': 17.34, 'i': 16.92}
errs = {'u': 0.03, 'g': 0.03, 'r': 0.03, 'i': 0.03}

mags_grid_file = h5py.File(work_dir+'/vst_mags_grid_z=%6.4f.hdf5'%redshift, 'r')

# prepares axes for interpolation
Z_grid = mags_grid_file['Z_grid'][()]
Z_spline = splrep(Z_grid, np.arange(len(Z_grid)), k=1, s=0)

tau_grid = mags_grid_file['tau_grid'][()]
tau_spline = splrep(tau_grid, np.arange(len(tau_grid)), k=3, s=0)

tau_V_grid = mags_grid_file['tau_V_grid'][()]
tau_V_spline = splrep(tau_V_grid, np.arange(len(tau_V_grid)), k=3, s=0)

age_grid = mags_grid_file['age_grid'][()]
age_spline = splrep(age_grid, np.arange(len(age_grid)), k=3, s=0)

axes = {0: Z_spline, 1: tau_spline, 2: tau_V_spline, 3: age_spline}

interp = {}
for band in bands:
    interp[band] = ndinterp.ndInterp(axes, mags_grid_file['%s_mag'%band][()])

# guesses stellar mass given reasonable values of other parameters
uniage = pygalaxev_cosmology.uniage(redshift)/1e9

start_pnt = {}
start_pnt['tau'] = 1.
start_pnt['age'] = 0.8 * uniage
start_pnt['Z'] = 0.02
start_pnt['tau_V'] = 0.1

par_list = ['Z', 'tau', 'tau_V', 'age']
start_arr = np.empty(4)
for i in range(4):
    start_arr[i] = start_pnt[par_list[i]]

mmags = []
mobs = []
merr = []
for band in bands:
    mmags.append(interp[band].eval(start_arr.reshape((1, 4)))[0])
    mobs.append(mags[band])
    merr.append(errs[band])
mmags = np.array(mmags)
mobs = np.array(mobs)
merr = np.array(merr)

# minimum chi^2 stellar mass
logmass_guess = 1/2.5 * ((mmags - mobs)/merr**2).sum() / (1./merr**2).sum()
print('starting value of log-stellar mass: %3.2f'%logmass_guess)

# defines allowed range of parameters and initial spread of walker positions
mass_par = {'name': 'logmass', 'lower': 8.5, 'upper': 13., 'guess': logmass_guess, 'spread': 0.1}
age_par = {'name': 'age', 'lower': 0., 'upper': uniage, 'guess': start_pnt['age'], 'spread': 1.}
tau_par = {'name': 'tau', 'lower': 0., 'upper': uniage, 'guess': start_pnt['tau'], 'spread': 1.}
logZ_par = {'name': 'logZ', 'lower': -4., 'upper': -1., 'guess': np.log10(start_pnt['Z']), 'spread': 0.1}
logtau_V_par = {'name': 'tau_V', 'lower': -2., 'upper': 0.3, 'guess': np.log10(start_pnt['tau_V']), 'spread': 0.1}

pars = [mass_par, age_par, tau_par, logZ_par, logtau_V_par]
npars = len(pars)

bounds = []
for par in pars:
    bounds.append((par['lower'], par['upper']))

def logprior(p):
    for i in range(npars):
        if p[i] < bounds[i][0] or p[i] > bounds[i][1]:
            return -1e300
    return 0.

def modelmags(p):
    mass, age, tau, logZ, logtau_V = p

    arr = np.array([10.**logZ, tau, 10.**logtau_V, age])

    mmags = []
    for band in bands:
        mmags.append(interp[band].eval(arr.reshape((1, 4)))[0] - 2.5*mass)

    return mmags 

def logpfunc(p):

    lprior = logprior(p)
    if lprior < 0.:
        return -1e300, -99.*np.ones(len(bands))

    mmags = modelmags(p)

    sumlogp = 0.

    for band, mmag in zip(bands, mmags):
        sumlogp += -0.5*(mags[band] - mmag)**2/errs[band]**2

    return sumlogp, mmags

sampler = emcee.EnsembleSampler(nwalkers, npars, logpfunc, threads=1)

start = []
for i in range(nwalkers):
    tmp = np.zeros(npars)
    for j in range(npars):
        a, b = (bounds[j][0] - pars[j]['guess'])/pars[j]['spread'], (bounds[j][1] - pars[j]['guess'])/pars[j]['spread']
        p0 = truncnorm.rvs(a, b, size=1)*pars[j]['spread'] + pars[j]['guess']
        tmp[j] = p0

    start.append(tmp)

print("Sampling")

sampler.run_mcmc(start, nstep, progress=True)

blobchain = sampler.blobs

output_file = h5py.File(outname, 'w')

saveind = thin*np.arange((nstep-burnin)//thin) + burnin

for n in range(npars):
    output_file.create_dataset(pars[n]['name'], data=sampler.chain[:, saveind, n])
output_file.create_dataset('logp', data=sampler.lnprobability[:, saveind])

magdic = {}
for band in bands:
    magdic['mag_%s'%band] = np.zeros((nwalkers, nstep))

for i in range(nstep):
    for j in range(nwalkers):
        for n in range(len(bands)):
            magdic['mag_%s'%bands[n]][j, i] = blobchain[i][j][n]

for band in bands:
    output_file.create_dataset('mag_%s'%band, data = magdic['mag_%s'%band][:, saveind])

