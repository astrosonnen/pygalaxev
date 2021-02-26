import numpy as np
from scipy.interpolate import splrep, splev


# combines filter transmission curve with
# - Atmospheric transmission
# - Mirror reflectivity
# - Primary focus unit throughput
# - Quantum efficiency of the detector
# to obtain transmission curves in the five HSC broad bands

make_plots = True # if True, makes .png files with the transmission curve of each filter

bands = ['g', 'r', 'i', 'z', 'y']

# reads in Quantum efficiency of CCD
f = open('hsc_response/qe_ccd_HSC.txt', 'r')
wav_qe, qe_arr = np.loadtxt(f, unpack=True)
f.close()

qe_spline = splrep(wav_qe, qe_arr, k=1) # spline interpolation

# reads in dewar window transmittance
f = open('hsc_response/throughput_win.txt', 'r')
wav_dew, dew_arr = np.loadtxt(f, unpack=True)
f.close()

dew_spline = splrep(wav_dew, dew_arr, k=1)

# primary focus unit
f = open('hsc_response/throughput_popt2.txt', 'r')
wav_popt2, popt2_arr = np.loadtxt(f, unpack=True)
f.close()

popt2_spline = splrep(wav_popt2, popt2_arr, k=1)

# reflectivity (using before 8th recoating on Oct 2017)
f = open('hsc_response/subaru_m1_r_20171010.txt', 'r')
wav_refl, refl_arr = np.loadtxt(f, usecols=(0, 2), unpack=True)
f.close()

refl_spline = splrep(np.flipud(wav_refl*10.), 0.01*np.flipud(refl_arr), k=1)

# atmospheric transmission
f = open('hsc_response/STD_BANDPASSES_DR1.dat.txt', 'r')
wav_atm, atm_arr = np.loadtxt(f, usecols=(0, 6), unpack=True)
f.close()

atm_spline = splrep(wav_atm, atm_arr, k=1)

for band in bands:
    f = open('hsc_response/HSC-%s.txt'%band, 'r')
    wav_filt, filt_arr = np.loadtxt(f, unpack=True)
    f.close()
    
    if make_plots:
        import pylab
        
        pylab.plot(wav_filt, filt_arr, label='Filter')
        pylab.plot(wav_filt, splev(wav_filt, qe_spline), label='QE')
        pylab.plot(wav_filt, splev(wav_filt, dew_spline), label='Dewar')
        pylab.plot(wav_filt, splev(wav_filt, popt2_spline), label='Primary focus unit')
        pylab.plot(wav_filt, splev(wav_filt, refl_spline), label='Mirror reflectivity')
        pylab.plot(wav_filt, splev(wav_filt, atm_spline), label='Atmosphere')
        
        full_transmission = filt_arr * splev(wav_filt, qe_spline) * splev(wav_filt, dew_spline) * splev(wav_filt, popt2_spline) * splev(wav_filt, refl_spline) * splev(wav_filt, atm_spline)
        full_transmission[full_transmission<0.] = 0.
        
        pylab.plot(wav_filt, full_transmission, label='Total')
        
        pylab.legend()
        pylab.ylim(-0.1, 1.1)
        
        pylab.savefig('HSC_%s.png'%band)
        pylab.close()
    
    f = open('HSC_%s.res'%band, 'w')
    np.savetxt(f, np.array((wav_filt, full_transmission)).T, fmt=('%2.1f', '%5.4f'))
    f.close()

