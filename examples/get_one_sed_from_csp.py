import pygalaxev
import numpy as np
from scipy.interpolate import splrep, splev
import os
import h5py


work_dir = '/data2/sonnenfeld/stellarpop_csp_models/bc03_basel_chabrier/'
grid_file = h5py.File(work_dir+'/BaSeL_Chabrier_sed_grid.hdf5', 'r')


