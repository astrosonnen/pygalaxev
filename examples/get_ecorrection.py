import os
import pygalaxev
import pygalaxev_cosmology
from pygalaxev_cosmology import c as csol, L_Sun, Mpc
import numpy as np
import h5py
from scipy.interpolate import splrep, splev, splint


# creates CSP models on a grid of stellar population parameters using galaxev

# selects the stellar template library:
# Low-resolution 'BaSeL' library, Chabrier IMF
ssp_dir = '/data2/sonnenfeld/bc03/BaSeL3.1_Atlas/Chabrier_IMF/'
tempname = 'lr_BaSeL'

work_dir = './'

# Using Padova 1994 tracks.
Z = 0.02
Zcode  ='m62'

mu = 0.3 # fraction of attenuation due to diffuse interstellar medium (fixed)
epsilon = 0. # gas recycling (no recycling if set to zero) (fixed)

nwav = 2023 # size of wavelength grid (can be looked up by running 'csp_galaxev' on any .ised file of the spectral library)

# fixed of values of tau and tau_V
tau = 1.
tau_V = 0.1

# formation redshift
z_form = 3.
age_z0 = pygalaxev_cosmology.lookback(z_form)/1e9
z_max = 1.
nz = 21
z_grid = np.linspace(0., z_max, nz)
z_grid = np.flipud(z_grid)
age_grid = 0.*z_grid
for i in range(nz):
    age_grid[i] = age_z0 - pygalaxev_cosmology.lookback(z_grid[i])/1e9
nage = nz

lum_grid = np.zeros(nz)

tmpname = work_dir+'/tmp.in'

# Create the models
isedname = ssp_dir+'/bc2003_%s_%s_chab_ssp.ised'%(tempname, Zcode)
cspname = 'bc03_chab_Z=%6.4f_tau=%5.3f_tV=%5.3f_mu=%3.1f_eps=%5.3f'%(Z, tau, tau_V, mu, epsilon)

pygalaxev.run_csp_galaxev(isedname, cspname, sfh_pars=tau, tau_V=tau_V, mu=0.3, epsilon=0., work_dir=work_dir)

# Create the mass normalization models
massname = work_dir+'/'+cspname+'.mass'
d = np.loadtxt(massname)
mass_spline = splrep(d[:, 0], d[:, 10], k=3, s=0) #using the sum of M*_liv+M_rem to renormalize the mass

# extracts SEDs on age grid
oname = work_dir+'/'+cspname+'_agegrid.sed'
pygalaxev.create_galaxevpl_config(tmpname, work_dir+'/'+cspname+'.ised', oname, age_grid)
os.system('$bc03/galaxevpl < %s'%tmpname)

f = open(oname, 'r')
wsed = np.loadtxt(f)
f.close()

wave = wsed[:, 0]

flux_grid = wsed[:, 1:]

# Renormalize the mass to 1M_Sun at z=0
logAge = np.log10(age_z0)+9.
mass = splev(logAge, mass_spline)
sed_grid = flux_grid/mass

Dlum = 1e-5 # luminosity distance in Mpc
M_r_Sun = 4.65 # r-band AB absolute magnitude of the Sun

# reads the filter response curve
filtdir = os.environ.get('PYGALAXEVDIR') + '/filters/'
filtname = filtdir+'r_SDSS.res'

f = open(filtname, 'r')
filt_wave, filt_t = np.loadtxt(f, unpack=True)
f.close()

filt_spline = splrep(filt_wave, filt_t)
wmin_filt, wmax_filt = filt_wave[0], filt_wave[-1]

cond_filt = (wave>=wmin_filt)&(wave<=wmax_filt)
nu_cond = np.flipud(cond_filt)

# Evaluate the filter response at the wavelengths of the spectrum
response = splev(wave[cond_filt], filt_spline)
nu_filter = csol*1e8/wave[cond_filt]

# flips arrays
response = np.flipud(response)
nu_filter = np.flipud(nu_filter)

# filter normalization
bp = splrep(nu_filter, response/nu_filter, s=0, k=1)
bandpass = splint(nu_filter[0], nu_filter[-1], bp)

# loops over the age grid, computes r-band Luminosity
for a in range(nage):
    llambda = sed_grid[:, a]

    flambda_obs = llambda*L_Sun/(4.*np.pi*(Dlum*Mpc)**2) # observed specific flux in erg/s/cm^2/AA
    fnu = flambda_obs * wave**2 / csol * 1e-8 # F_nu in cgs units

    nu_obs = np.flipud(csol/wave*1e8)
    fnu = np.flipud(fnu)

    # Integrate
    observed = splrep(nu_filter, response*fnu[nu_cond]/nu_filter, s=0, k=1)
    flux = splint(nu_filter[0], nu_filter[-1], observed)

    M_r = -2.5*np.log10(flux/bandpass) -48.6

    lum_grid[a] = 10.**(-(M_r - M_r_Sun)/2.5)

# Clean up
os.system('rm %s'%oname)

grid_file = h5py.File('ecorr_rband_Chab_zform3.0_Z=0.02.hdf5', 'w')
grid_file.attrs['units'] = 'r-band L_Sun for 1M_Sun (living + remnants) at z=0'
grid_file.attrs['tau'] = tau
grid_file.attrs['tau_V'] = tau_V
grid_file.attrs['Z'] = Z
grid_file.attrs['z_form'] = z_form

grid_file.create_dataset('age_grid', data=age_grid)
grid_file.create_dataset('redshift_grid', data=z_grid)
grid_file.create_dataset('lum_grid', data=lum_grid)

