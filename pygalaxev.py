import os
import numpy as np


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

    if sfh == 'tau':
        tau = sfh_pars
    else:
        raise ValueError("Only 'tau' SFH is currently implemented")

    f = open(tmpname, 'w')

    # Metallicity/IMF
    f.write('%s\n'%isedname)
    # Use dust with tau_V = tV and mu = mu
    f.write('Y\n')
    f.write('%f\n'%tau_V)
    f.write('%f\n'%mu)
    f.write('0\n') # don't compute flux-weighted age
    # Use exponential SFH with tau = t
    f.write('1\n')
    f.write('%f\n'%tau)
    # Choose whether or not to recycle gas
    if epsilon!=0:
        f.write('Y\n')
        f.write('%f\n'%epsilon)
    else:
        f.write('N\n')
    # Cutoff star formation at 20 Gyr
    f.write('20\n')
    f.write('%s/%s\n'%(work_dir, output_tmpname))
    f.close()

    # Run bc03
    os.system('$bc03/csp_galaxev < %s'%tmpname)

    # Keep the output *.ised and *.4color files
    os.system('cp %s/%s.ised %s/%s.ised'%(work_dir, output_tmpname, work_dir, outname))
    os.system('cp %s/%s.4color %s/%s.mass'%(work_dir, output_tmpname, work_dir, outname))

    # Clean up
    os.system('rm -f %s/%s*'%(work_dir, output_tmpname))

