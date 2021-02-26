import numpy as np
from scipy.interpolate import splrep, splev
import pylab


filtdir = 'VST_response/'

# optical atmospheric transmission 
f = open('paranal_atm_transmission.dat', 'r')
atm_wav, atm_t = np.loadtxt(f, unpack=True)
f.close()

atm_wav *= 10.

atm_spline = splrep(atm_wav, atm_t, k=1)

pylab.plot(atm_wav, atm_t, color='k', label='Atmosphere')

bands = ['u', 'g', 'r', 'i']

colors = ['b', 'c', 'g', 'orange', 'r']

for band, color in zip(bands, colors):

    f = open('%s_OmegaCAM.res'%band, 'r')
    old_wav, old_t = np.loadtxt(f, unpack=True)
    f.close()

    f = open(filtdir+'/%s_modified_transmission_data_with_ccd'%band, 'r')
    filt_tab = np.loadtxt(f)
    f.close()

    filt_wav = filt_tab[:, 0]
    filt_t = filt_tab[:, 1]

    filt_wav *= 10.

    cut = (filt_wav >= old_wav[0]) & (filt_wav <= old_wav[-1])

    filt_wav = filt_wav[cut]
    filt_t = filt_t[cut]

    full_transmission = filt_t * splev(filt_wav, atm_spline)
    full_transmission[full_transmission<0.] = 0.
        
    pylab.plot(filt_wav, full_transmission, label=band, color=color)
   
    f = open('OmegaCAM_watm_%s.res'%band, 'w')
    np.savetxt(f, np.array((filt_wav, full_transmission)).T, fmt=('%2.1f', '%5.4f'))
    f.close()

pylab.legend()
pylab.xlim(3000., 10000.)
pylab.ylim(-0.1, 1.1)
        
pylab.savefig('OmegaCAM_filters.png')
pylab.show()
 
