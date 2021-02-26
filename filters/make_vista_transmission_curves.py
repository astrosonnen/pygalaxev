import numpy as np
from scipy.interpolate import splrep, splev
import pylab


filtdir = 'VISTA_Filters_QE_Atm_curves/'

f = open(filtdir+'trans_10_10.dat', 'r')
atm1_wav, atm1_t = np.loadtxt(f, unpack=True)
f.close()

atm1_wav *= 1e4

f = open(filtdir+'trans_30_20.dat', 'r')
atm2_wav, atm2_t = np.loadtxt(f, unpack=True)
f.close()

atm2_wav *= 1e4

# optical atmospheric transmission (from HSC, because the VISTA one doesn't go blue enough)
f = open('hsc_response/STD_BANDPASSES_DR1.dat.txt', 'r')
opt_atm_wav, opt_atm_t = np.loadtxt(f, usecols=(0, 6), unpack=True)
f.close()

ir_atm_spline = splrep(atm1_wav, atm1_t, k=1)
opt_atm_spline = splrep(opt_atm_wav, opt_atm_t, k=1)
atm_splines = {'optical': opt_atm_spline, 'ir': ir_atm_spline}

pylab.plot(atm1_wav, atm1_t, color='k', label='Atmosphere (IR)')
#pylab.plot(atm2_wav, atm2_t, color='k')
pylab.plot(opt_atm_wav, opt_atm_t, color='grey', label='Atmosphere (optical)')

# Mirror reflectivities
f = open(filtdir+'/VISTA_M1_Reflectivity_forETC_2009-Sep12.txt', 'r')
m1_wav, m1_t = np.loadtxt(f, unpack=True)
f.close()

m1_wav *= 10.
m1_t /= 100.

mirror_spline = splrep(m1_wav, m1_t)

f = open(filtdir+'/VISTA_M2_Reflectivity_forETC_2007-Jun19.txt', 'r')
m2_wav, m2_t = np.loadtxt(f, unpack=True)
f.close()

m2_wav *= 10.
m2_t /= 100.

pylab.plot(m1_wav, m1_t, color='c', label='Mirror reflectivity')
#pylab.plot(m2_wav, m2_t, color='grey')

# Quantum efficiency
f = open(filtdir+'/qe.tab', 'r')
qe_wav, qe_t = np.loadtxt(f, unpack=True)
f.close()

qe_wav *= 10.
qe_t /= 100.

qe_spline = splrep(qe_wav, qe_t, k=1)

pylab.plot(qe_wav, qe_t, color='orange', label='Quantum efficiency')

bands = ['Z', 'Y', 'J', 'H', 'Ks']
regions = [(qe_wav[0], 10000.), (9000., 11500.), (11000., 14500.), (14000., 19000.), (18500., 24500.)]
colors = ['b', 'g', 'y', 'r', 'brown']

for band, reg, color in zip(bands, regions, colors):

    f = open(filtdir+'/VISTA_Filters_at80K_forETC_%s.dat'%band, 'r')
    filt_wav, filt_t = np.loadtxt(f, unpack=True)
    f.close()

    filt_wav *= 10.
    filt_t /= 100.

    cut = (filt_wav > reg[0]) & (filt_wav < reg[1])

    filt_wav = filt_wav[cut]
    filt_t = filt_t[cut]

    if band == 'Z':
        atm_transmission = splev(filt_wav, atm_splines['optical'])
    else:
        atm_transmission = splev(filt_wav, atm_splines['ir'])

    full_transmission = filt_t * splev(filt_wav, qe_spline) * splev(filt_wav, mirror_spline) * atm_transmission
    full_transmission[full_transmission<0.] = 0.
        
    pylab.plot(filt_wav, full_transmission, label=band, color=color)
   
    f = open('VISTA_%s.res'%band, 'w')
    np.savetxt(f, np.array((filt_wav, full_transmission)).T, fmt=('%2.1f', '%5.4f'))
    f.close()

pylab.legend()
pylab.xlim(8000., 25000.)
pylab.ylim(-0.1, 1.1)
        
pylab.savefig('VISTA_filters.png')
pylab.close()
 
