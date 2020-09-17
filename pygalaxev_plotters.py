#series of plotting tools

import pylab
import numpy as np
#from matplotlib import rc
#rc('text',usetex=True)

def rgb_alpha(color,alpha):
    if type(color[0]) == int:
        fcolor = []
        for col in color:
            fcolor.append(col/255.)
    else:
        fcolor = color
        
    R = alpha*fcolor[0] + 1. - alpha
    G = alpha*fcolor[1] + 1. - alpha
    B = alpha*fcolor[2] + 1. - alpha

    return (R,G,B)

def rgb_to_hex(color):
    digits = ['0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f']
    hexx = '#'
    for col in color:
        col = int(col*255)
        first = digits[col//16]
        second = digits[col%16]
        hexx += first
        hexx += second
    return hexx

def probcontour(xarr,yarr,style='lines',color='k',smooth=2,bins=100,weights=None, linewidths=1, linestyles='solid', nlevels=3):
    from scipy import ndimage
    import numpy
    H,xbins,ybins = pylab.histogram2d(xarr,yarr,bins=bins,weights=weights)

    H = ndimage.gaussian_filter(H,smooth)
    sortH = numpy.sort(H.flatten())
    cumH = sortH.cumsum()
# 1, 2, 3-sigma, for the old school:
    lvl68 = sortH[cumH>cumH.max()*0.32].min()
    lvl95 = sortH[cumH>cumH.max()*0.05].min()
    lvl99 = sortH[cumH>cumH.max()*0.003].min()
    if style == 'lines':
        #pylab.contour(H.T,[lvl99,lvl95,lvl68],colors=color,\
        #                  extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]), linewidths=linewidths, linestyles=linestyles)
        pylab.contour(H.T,[lvl95,lvl68],colors=color,\
                          extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]), linewidths=linewidths, linestyles=linestyles)


    else:
        # Shaded areas:
        if type(color) == str:
            pylab.contourf(H.T,[lvl95,lvl68],colors=color,alpha=0.5,\
                               extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]))
            pylab.contourf(H.T,[lvl68,1e8],colors=color,alpha=0.9,\
                               extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]))
        else:
            pylab.contourf(H.T,[lvl95,lvl68],colors=rgb_to_hex(rgb_alpha(color,0.5)),\
                               extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]))
            pylab.contourf(H.T,[lvl68,1e8],colors=rgb_to_hex(rgb_alpha(color,0.9)),\
                               extent=(xbins[0],xbins[-1],ybins[0],ybins[-1]))

