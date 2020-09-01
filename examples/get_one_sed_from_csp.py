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

age = 13. # time since beginning of star formation (in Gyr)
wrange = None # wavelength range (if None, uses whole range provided by library)

cspname = 'bc03_Z=%6.4f_tau=%5.3f_tV=%5.3f_mu=%3.1f_eps=%5.3f'%(Z, tau, tau_V, mu, epsilon)
tmpname = 'tmp.in'
oname = cspname + '_age=%06.3f.sed'%age

pygalaxev.create_galaxevpl_config(work_dir+'/'+tmpname, work_dir+'/'+cspname+'.ised', work_dir+'/'+oname, age, wrange)
os.system('$bc03/galaxevpl < %s'%tmpname)

f = open(work_dir+'/'+oname, 'r')
wsed = np.loadtxt(f)
f.close()

wave = wsed[:, 0]
flux = wsed[:, 1]

# Create the mass normalization models
massname = cspname+'.mass'
d = np.loadtxt(massname)
mass_spline = splrep(d[:, 0], d[:, 10], k=3, s=0) #using the sum of M*_liv+M_rem to renormalize the mass

# Renormalize the mass!
logAge = np.log10(age)+9.
mass = splev(logAge, mass_spline)
sed = flux/mass

# Clean up
os.system('rm %s/%s'%(work_dir, oname))

# writes SED to a separate file
sedfile = h5py.File(work_dir+'/'+oname.replace('.sed', '.hdf5'), 'w')

sedfile.create_dataset('wave', data=wave)
dset = sedfile.create_dataset('Llambda', data=flux)
dset.attrs['units'] = 'Llambda (in units of L_Sun/Angstrom) for 1M_Sun (living + remnants)'

sedfile.close()

