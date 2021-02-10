import numpy as np
import h5py
import pygalaxev_cosmology
from pygalaxev_cosmology import L_Sun, Mpc, c as csol
from scipy.interpolate import splrep, splev, splint
import os


# this code reads in a CSP SED file and, given a galaxy redshift, calculates the magnitudes in a series of filters

work_dir = './'
pygalaxevdir = os.environ.get('PYGALAXEVDIR')
filtdir = pygalaxevdir+'/filters/'

redshift = 0.1945
Dlum = pygalaxev_cosmology.Dlum(redshift) # luminosity distance in Mpc

log_mstar = 11.

Z = 0.02
tau = 1.
tau_V = 0.1
mu = 0.3
epsilon = 0.

age = 11. # time since beginning of star formation (in Gyr)

cspname = 'bc03_Z=%6.4f_tau=%5.3f_tV=%5.3f_mu=%3.1f_eps=%5.3f'%(Z, tau, tau_V, mu, epsilon)
sedname = cspname + '_age=%06.3f.hdf5'%age

sedfile = h5py.File(work_dir+'/'+sedname, 'r')

wave = sedfile['wave'][()]
llambda = sedfile['Llambda'][()]

bands = ['u', 'g', 'r', 'i']

wave_obs = wave * (1.+redshift)
flambda_obs = llambda*L_Sun/(4.*np.pi*(Dlum*Mpc)**2)/(1.+redshift) # observed specific flux in erg/s/cm^2/AA
fnu = flambda_obs * wave_obs**2 / csol * 1e-8 # F_nu in cgs units

nu_obs = np.flipud(csol/wave_obs*1e8)
fnu = np.flipud(fnu)

for band in bands:
    # loads filter transmission curve file
    filtname = filtdir+'/%s_OmegaCAM.res'%band

    f = open(filtname, 'r')
    filt_wave, filt_t = np.loadtxt(f, unpack=True)
    f.close()

    filt_spline = splrep(filt_wave, filt_t)

    wmin_filt, wmax_filt = filt_wave[0], filt_wave[-1]
    cond_filt = (wave_obs>=wmin_filt)&(wave_obs<=wmax_filt)
    nu_cond = np.flipud(cond_filt)

    # Evaluate the filter response at the wavelengths of the spectrum
    response = splev(wave_obs[cond_filt], filt_spline)
    nu_filter = csol*1e8/wave_obs[cond_filt]

    # flips arrays
    response = np.flipud(response)
    nu_filter = np.flipud(nu_filter)

    # filter normalization
    bp = splrep(nu_filter, response/nu_filter, s=0, k=1)
    bandpass = splint(nu_filter[0], nu_filter[-1], bp)

    # Integrate
    observed = splrep(nu_filter, response*fnu[nu_cond]/nu_filter, s=0, k=1)
    flux = splint(nu_filter[0], nu_filter[-1], observed)

    mag = -2.5*np.log10(flux/bandpass) -48.6 -2.5*log_mstar
    print(mag)
 
