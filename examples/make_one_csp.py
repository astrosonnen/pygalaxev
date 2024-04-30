import os
import pygalaxev


# creates one CSP model using galaxev

# selects the stellar template library:

# Low-resolution 'BaSeL' library, Chabrier IMF
ssp_dir = '/Users/alessandro/science/bc03/BaSeL3.1_Atlas/Chabrier_IMF/'
tempname = 'lr_BaSeL'
work_dir = './'

# Using Padova 1994 tracks. Selects the metallicity:

Z = 0.02 # solar metallicity
Zcode_dic = {0.0001:'m22', 0.0004:'m32', 0.004:'m42', 0.008:'m52', 0.02:'m62', 0.05:'m72', 0.1:'m82'}
Zcode = Zcode_dic[Z]

tau = 1. # exponential star formation history timescale (in Gyr)
tau_V = 0.1 # dust effective optical depth
mu = 0.3 # fraction of attenuation due to diffuse interstellar medium
epsilon = 0. # gas recycling (no recycling if set to zero)

isedname = ssp_dir+'/bc2003_%s_%s_chab_ssp.ised'%(tempname, Zcode)
outname = 'bc03_Z=%6.4f_tau=%5.3f_tV=%5.3f_mu=%3.1f_eps=%5.3f'%(Z, tau, tau_V, mu, epsilon)

pygalaxev.run_csp_galaxev(isedname, outname, sfh_pars=tau, mu=0.3, epsilon=0., work_dir=work_dir)

