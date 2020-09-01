import numpy as np


filtdir = '/data1/sonnenfeld/downloads/Filters_QE_Atm_curves/'

bands = ['Z']

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

