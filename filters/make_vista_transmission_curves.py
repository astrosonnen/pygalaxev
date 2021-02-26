import numpy as np
import pylab
from scipy.interpolate import splrep, splev


filtdir = 'VISTA_Filters_QE_Atm_curves/'

f = open(filtdir+'trans_10_10.dat', 'r')
atm1_wav, atm1_t = np.loadtxt(f, unpack=True)
f.close()

f = open(filtdir+'trans_30_20.dat', 'r')
atm2_wav, atm2_t = np.loadtxt(f, unpack=True)
f.close()

pylab.plot(atm1_wav, atm1_t)
pylab.plot(atm2_wav, atm2_t)
pylab.show()

bands = ['Z']

"""
for band in bands:

    f = open(filtdir+'/VISTA_Filters_at80K_forETC_%s.dat'%band, 'r')
    wav, t = np.loadtxt(f, unpack=True)
    f.close()

    t[t<0.] = 0.
    t /= 100.

    cut = (wav > 700.) & (wav < 1000.)

    f = open('VISTA_%s.res'%band, 'w')
    np.savetxt(f, np.array((10.*wav[cut], t[cut])).T, fmt='%f')
    f.close()
"""

