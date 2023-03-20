import os
import numpy as np
import pygalaxev_cosmology
from pygalaxev_cosmology import Mpc, c as csol, L_Sun
from scipy.interpolate import splrep, splev, splint


def create_galaxevpl_config(configname, cspname, outname, age, wrange=None):
    age = np.atleast_1d(age)
    age_str = '%4.3f'%age[0]
    if len(age) > 1:
        for i in range(len(age)-1):
            age_str += ', %4.3f'%age[i+1]

    if wrange is not None:
        wrange_str = '%f, %f'%(wrange[0], wrange[1])
    else:
        wrange_str = ''

    f = open(configname, 'w')
    f.write('%s\n'%cspname)
    f.write('%s\n'%age_str)
    f.write('%s\n'%wrange_str)
    f.write('%s\n'%outname)
    f.close()

def run_csp_galaxev(isedname, outname, sfh='tau', sfh_pars=1., tau_V=0.1, mu=0.3, epsilon=0., work_dir='./', input_tmpname='tmp.in', output_tmpname='mySSP'):
    """
    sfh: star formation history model. Exponentially declining, by default.
    sfh_pars: parameter describing the star formation history. If sfh='tau', then sfh_pars is the timescale of the exponential (in Gyr), indicated as tau.
    tau_V: dust optical depth
    mu: fraction of dust due to diffuse interstellar medium
    epsilon: gas recycling
    """

    tmpname = work_dir+'/%s'%input_tmpname

    # prepares input file for csp_galaxev
    inputlines = []
    # Metallicity/IMF
    inputlines.append('%s\n'%isedname)
    # Use dust with tau_V = tV and mu = mu
    inputlines.append('Y\n')
    inputlines.append('%f\n'%tau_V)
    inputlines.append('%f\n'%mu)
    inputlines.append('0\n') # don't compute flux-weighted age
    # choose star formation history
    if sfh == 'tau':
        inputlines.append('1\n')
        tau = sfh_pars
        inputlines.append('%f\n'%tau)
        # Choose whether or not to recycle gas
        if epsilon!=0:
            inputlines.append('Y\n')
            inputlines.append('%f\n'%epsilon)
        else:
            inputlines.append('N\n')
        # Cutoff star formation at 20 Gyr
        inputlines.append('20\n')
    elif sfh == 'SSP':
        inputlines.append('0\n')
    else:
        raise ValueError("Only 'tau' and 'SSP' SFHs are currently implemented")

    inputlines.append('%s/%s\n'%(work_dir, output_tmpname))

    f = open(tmpname, 'w')
    f.writelines(inputlines)
    f.close()

    # Run bc03
    os.system('$bc03/csp_galaxev < %s'%tmpname)

    # Keep the output *.ised and *.4color files
    os.system('cp %s/%s.ised %s/%s.ised'%(work_dir, output_tmpname, work_dir, outname))
    os.system('cp %s/%s.4color %s/%s.mass'%(work_dir, output_tmpname, work_dir, outname))

    # Clean up
    os.system('rm -f %s/%s*'%(work_dir, output_tmpname))

def get_mag_from_sed(wave, llambda, redshift, filtname, cosmo=pygalaxev_cosmology.default_cosmo):

    filtdir = os.environ.get('PYGALAXEVDIR') + '/filters/'
    if redshift == 0.: # computes absolute magnitude if z=0
        Dlum = 1e-5
    else:
        Dlum = pygalaxev_cosmology.Dlum(redshift, cosmo=cosmo) # luminosity distance in Mpc

    wave_obs = wave * (1.+redshift)
    flambda_obs = llambda*L_Sun/(4.*np.pi*(Dlum*Mpc)**2)/(1.+redshift) # observed specific flux in erg/s/cm^2/AA
    fnu = flambda_obs * wave_obs**2 / csol * 1e-8 # F_nu in cgs units

    nu_obs = np.flipud(csol/wave_obs*1e8)
    fnu = np.flipud(fnu)

    fullfiltname = filtdir+filtname

    f = open(fullfiltname, 'r')
    filt_wave, filt_t = np.loadtxt(f, unpack=True)
    f.close()

    filt_spline = splrep(filt_wave, filt_t)

    wmin_filt, wmax_filt = filt_wave[0], filt_wave[-1]
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

    # Integrate
    observed = splrep(nu_filter, response*fnu[nu_cond]/nu_filter, s=0, k=1)
    flux = splint(nu_filter[0], nu_filter[-1], observed)

    mag = -2.5*np.log10(flux/bandpass) -48.6

    return mag

