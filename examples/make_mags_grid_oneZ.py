import os
import h5py
import pygalaxev
import pygalaxev_cosmology
from pygalaxev_cosmology import c as csol, L_Sun, Mpc
from scipy.interpolate import splrep, splev, splint
import numpy as np


# calculates VST model magnitudes on the grid of stellar population parameters created with the 'make_sed_grid.py' script

work_dir = '/data2/sonnenfeld/stellarpop_csp_models/bc03_basel_chabrier/'
grid_file = h5py.File(work_dir+'/BaSeL_Chabrier_sed_grid.hdf5', 'r')
filtdir = os.environ.get('PYGALAXEVDIR') + '/filters/'

bands = ['u', 'g', 'r', 'i']

Z_fixed = 0.02
Zgrid_ind = 4
output_file = h5py.File(work_dir+'/vst_mags_grid_Z=%3.2f_zgrid.hdf5'%Z_fixed, 'w')

wave = grid_file['wave'][()]

sed_grid = grid_file['sed_grid'][()]
nZ, ntau, ntau_V, nage, nwave = sed_grid.shape

dz = 0.01
redshift_grid = np.arange(0.01, 1.00, dz)
nredshift = len(redshift_grid)

for band in bands:
    filtname = filtdir+'%s_OmegaCAM.res'%band

    f = open(filtname, 'r')
    filt_wave, filt_t = np.loadtxt(f, unpack=True)
    f.close()

    filt_spline = splrep(filt_wave, filt_t)
    wmin_filt, wmax_filt = filt_wave[0], filt_wave[-1]

    mag_grid = np.zeros((nredshift, ntau, ntau_V, nage))

    for i in range(nredshift):
        Dlum = pygalaxev_cosmology.Dlum(redshift_grid[i]) # luminosity distance in Mpc

        wave_obs = wave * (1.+redshift_grid[i])

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

        print('%s band, redshift %3.2f'%(band, redshift_grid[i]))
        for j in range(ntau):
            for k in range(ntau_V):
                for l in range(nage):
                    llambda = sed_grid[Zgrid_ind, j, k, l, :]

                    flambda_obs = llambda*L_Sun/(4.*np.pi*(Dlum*Mpc)**2)/(1.+redshift_grid[i]) # observed specific flux in erg/s/cm^2/AA
                    fnu = flambda_obs * wave_obs**2 / csol * 1e-8 # F_nu in cgs units

                    nu_obs = np.flipud(csol/wave_obs*1e8)
                    fnu = np.flipud(fnu)

                    # Integrate
                    observed = splrep(nu_filter, response*fnu[nu_cond]/nu_filter, s=0, k=1)
                    flux = splint(nu_filter[0], nu_filter[-1], observed)

                    mag_grid[i, j, k, l] = -2.5*np.log10(flux/bandpass) -48.6

    mag_dset = output_file.create_dataset('%s_mag'%band, data=mag_grid)

    mag_dset.attrs['axis_0'] = 'redshift'
    mag_dset.attrs['axis_1'] = 'tau'
    mag_dset.attrs['axis_2'] = 'tau_V'
    mag_dset.attrs['axis_3'] = 'age'
    mag_dset.attrs['axis_4'] = 'Wavelength'

output_file.create_dataset('redshift_grid', data=redshift_grid)
output_file.create_dataset('tau_grid', data=grid_file['tau_grid'][()])
output_file.create_dataset('tau_V_grid', data=grid_file['tau_V_grid'][()])
output_file.create_dataset('age_grid', data=grid_file['age_grid'][()])

