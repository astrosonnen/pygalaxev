import pylab
import numpy as np
import h5py
from pygalaxev_plotters import probcontour
import os
from matplotlib.ticker import MultipleLocator
from matplotlib import rc
rc('text', usetex=True)


fsize = 12

color = 'r'

chainname = 'mcmc_chain.hdf5'
chain_file = h5py.File(chainname, 'r')

print(chain_file.keys())

nbins = 20
burnin = 0

smooth = 1

fsize = 12

pars = ['logmass', 'logZ', 'tau', 'age', 'tau_V']

chain = {}
for par in pars:
    chain[par] = chain_file[par][:, burnin:].flatten()

labels = ['$\log{M_*}$', '$\log{Z}$', '$\\tau$ (Gyr)', 'age (Gyr)', '$\log{\\tau_V}$']

lims = [(10.9, 11.7), (-1.75, -1.35), (0., 2.), (0., 9.), (-1.6, 0.3)]
major_step = [0.5, 0.2, 1, 5, 1]
minor_step = [0.1, 0.05, 0.5, 1, 0.2]

npars = len(pars)

fig = pylab.figure()
pylab.subplots_adjust(left=0.1, right=0.99, bottom=0.1, top=0.99, hspace=0.1, wspace=0.1)

for i in range(npars):

    ax = fig.add_subplot(npars, npars, (npars+1)*i + 1)

    bins = np.linspace(lims[i][0], lims[i][1], nbins+1)

    pylab.hist(chain[pars[i]], bins=bins, color=color)

    if i==0:
        ylim = pylab.ylim()
        pylab.ylim(ylim[0], ylim[1])

    ax.set_xlim((lims[i][0], lims[i][1]))
    ax.tick_params(which='both', direction='in')
    ax.xaxis.set_major_locator(MultipleLocator(major_step[i]))
    ax.xaxis.set_minor_locator(MultipleLocator(minor_step[i]))

    ax.set_yticks(())
    if i == npars-1:
        ax.set_xlabel(labels[i], fontsize=fsize)
    else:
        ax.tick_params(axis='x', labelbottom=False)

donelabel = False
for j in range(1, npars): # loops over rows
    if j == npars-1:
        xvisible = True
    else:
        xvisible = False

    for i in range(j): # loops over columns
        ax = pylab.subplot(npars, npars, npars*j+i+1)

        xbins = np.linspace(lims[i][0], lims[i][1], nbins+1)
        ybins = np.linspace(lims[j][0], lims[j][1], nbins+1)

        probcontour(chain[pars[i]], chain[pars[j]], bins=(xbins, ybins), smooth=smooth, style='filled', color=color, linewidths=2)

        ax.set_xlim(lims[i])
        ax.set_ylim(lims[j])
        ax.tick_params(which='both', direction='in', labelsize=fsize)

        ax.xaxis.set_major_locator(MultipleLocator(major_step[i]))
        ax.xaxis.set_minor_locator(MultipleLocator(minor_step[i]))

        ax.yaxis.set_major_locator(MultipleLocator(major_step[j]))
        ax.yaxis.set_minor_locator(MultipleLocator(minor_step[j]))

        if i == 0:
            yvisible = True
            ax.set_ylabel(labels[j], fontsize=fsize)

        else:
            ax.tick_params(axis='y', labelleft=False)

        if xvisible:
            ax.set_xlabel(labels[i], fontsize=fsize)
        else:
            ax.tick_params(axis='x', labelbottom=False)

pylab.show()

