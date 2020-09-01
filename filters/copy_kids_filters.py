import numpy as np


filtdir = '/disks/shear10/mj/dr4_sdss/filters/'

bands = ['u', 'g', 'r', 'i']

for band in bands:
    f = open(filtdir+'/%s.dat'%band, 'r')
    lines = f.readlines()
    f.close()

    newlines = []
    for line in lines:
        wav, t = line.split()
        newlines.append('%4.2f %f\n'%(10.*float(wav), float(t)))

    f = open('OmegaCAM_%s.res'%band, 'w')
    f.writelines(newlines)
    f.close()

