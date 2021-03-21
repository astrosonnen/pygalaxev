from scipy.interpolate import splrep, splev
import numpy as np
import os


Z_Sun = 0.02 # solar metallicity

thisdir = os.path.dirname(os.path.abspath(__file__))

f = open(thisdir+'/gallazzi05_table2.txt', 'r')
gall_table = np.loadtxt(f)
f.close()

mbins = gall_table[:, 0]
z16 = gall_table[:, 2]
z84 = gall_table[:, 3]

zcen = 0.5 * (z16 + z84)
zsig = 0.5 * (z84 - z16)

zcen_spline = splrep(mbins, zcen, k=3)
zsig_spline = splrep(mbins, zsig, k=3)

def mzprior(mchab):

    return (splev(mchab, zcen_spline) + np.log10(Z_Sun), splev(mchab, zsig_spline))

