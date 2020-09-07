import pygalaxev
import numpy as np
from scipy.interpolate import splrep, splev
import os
import h5py


work_dir = './'

Z = 0.02
tau = 1.
tau_V = 0.1
mu = 0.3
epsilon = 0.

cspname = 'bc03_Z=%6.4f_tau=%5.3f_tV=%5.3f_mu=%3.1f_eps=%5.3f'%(Z, tau, tau_V, mu, epsilon)

age = 13. # time since beginning of star formation (in Gyr)

# Create the mass normalization models
massname = work_dir+'/'+cspname+'.mass'
d = np.loadtxt(massname)
mass_spline = splrep(d[:, 0], d[:, 10], k=3, s=0) #using the sum of M*_liv+M_rem to renormalize the mass

# extracts SED corresponding to the given age

tmpname = work_dir+'/tmp.in'

oname = work_dir+'/'+cspname+'_age=%06.3f.sed'%age

pygalaxev.create_galaxevpl_config(tmpname, work_dir+'/'+cspname+'.ised', oname, age)
os.system('$bc03/galaxevpl < %s'%tmpname)

f = open(oname, 'r')
wsed = np.loadtxt(f)
f.close()

wave = wsed[:, 0]
flux = wsed[:, 1]
           
# Renormalize the mass!
logAge = np.log10(age)+9.
mass = splev(logAge, mass_spline)
sed = flux/mass

# store SED to an .hdf5 file

output_file = h5py.File(oname.replace('.sed', '.hdf5'), 'w')
wave_dset = output_file.create_dataset('wave', data=wave)
wave_dset.attrs['units'] = 'Angstrom'

sed_dset = output_file.create_dataset('llambda', data=wave)
sed_dset.attrs['description'] = 'Luminosity density (dL/dlambda) for a 1 Solar Mass (living + remnants) stellar population'
sed_dset.attrs['units'] = 'L_Sun/Angstrom'

