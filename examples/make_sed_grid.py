import os
import pygalaxev
import numpy as np
import h5py
from scipy.interpolate import splrep, splev


# creates CSP models on a grid of stellar population parameters using galaxev

# selects the stellar template library:
# Low-resolution 'BaSeL' library, Chabrier IMF
ssp_dir = '/data2/sonnenfeld/bc03/BaSeL3.1_Atlas/Chabrier_IMF/'
tempname = 'lr_BaSeL'

work_dir = '/data2/sonnenfeld/stellarpop_csp_models/bc03_basel_chabrier/'

# Using Padova 1994 tracks.
Z_given_code = {'m22':0.0001, 'm32':0.0004, 'm42':0.004, 'm52':0.008, 'm62': 0.02, 'm72': 0.05, 'm82': 0.1}

mu = 0.3 # fraction of attenuation due to diffuse interstellar medium (fixed)
epsilon = 0. # gas recycling (no recycling if set to zero) (fixed)

nwav = 2023 # size of wavelength grid (can be looked up by running 'csp_galaxev' on any .ised file of the spectral library)

nZ = len(Z_given_code) # size of metallicity grid
ntau = 35 # size of grid in exponential timescale
ntau_V = 21 # size of grid in dust attenuation
nage = 31 # size of grid in age

# grid of values of tau, tau_V and age
tau_grid = np.linspace(0.04, 9.1, ntau)
tau_V_grid = np.logspace(-2., np.log10(2.), ntau_V)
age_grid = np.linspace(0.6, 13.5, nage)
Z_grid = []

grid = np.zeros((nZ, ntau, ntau_V, nage, nwav))

tmpname = work_dir+'/tmp.in'

for m in range(2, 9): # loop over metallicities
    Zcode = 'm%d2'%m
    Z = Z_given_code[Zcode]
    Z_grid.append(Z)

    for t in range(ntau):
        for tV in range(ntau_V):
            # Create the models
            isedname = ssp_dir+'/bc2003_%s_%s_chab_ssp.ised'%(tempname, Zcode)
            cspname = 'bc03_Z=%6.4f_tau=%5.3f_tV=%5.3f_mu=%3.1f_eps=%5.3f'%(Z, tau_grid[t], tau_V_grid[tV], mu, epsilon)

            pygalaxev.run_csp_galaxev(isedname, cspname, sfh_pars=tau_grid[t], tau_V=tau_V_grid[tV], mu=0.3, epsilon=0., work_dir=work_dir)

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

            for a in range(nage):
                flux = wsed[:, a+1]
           
                # Renormalize the mass!
                logAge = np.log10(age_grid[a])+9.
                mass = splev(logAge, mass_spline)
                sed = flux/mass

                grid[m-2, t, tV, a, :] = sed

            # Clean up
            os.system('rm %s'%oname)

Z_grid = np.array(Z_grid)

grid_file = h5py.File(work_dir+'/BaSeL_Chabrier_sed_grid.hdf5', 'w')
grid_dset = grid_file.create_dataset('sed_grid', data=grid)
grid_dset.attrs['units'] = 'Llambda (in units of L_Sun/Angstrom) for 1M_Sun (living + remnants)'
grid_dset.attrs['axis_0'] = 'Metallicity'
grid_dset.attrs['axis_1'] = 'tau'
grid_dset.attrs['axis_2'] = 'tau_V'
grid_dset.attrs['axis_3'] = 'age'
grid_dset.attrs['axis_4'] = 'Wavelength'

grid_file.create_dataset('Z_grid', data=Z_grid)
grid_file.create_dataset('tau_grid', data=tau_grid)
grid_file.create_dataset('tau_V_grid', data=tau_V_grid)
grid_file.create_dataset('age_grid', data=age_grid)
grid_file.create_dataset('wave', data=wave)

