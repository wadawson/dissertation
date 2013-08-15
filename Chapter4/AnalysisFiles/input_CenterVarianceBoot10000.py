'''
This is an input script for running CenterVarianceBoot.py
'''
from CenterVarianceBoot import *

# Galaxy Density Input
# User inputs:
id_sci = 0
N_boot = 10000
id_boot = 3

r_step = 1
r_stop = 21 #~200 kpc
N_bins = 50 #number of bins or histogram plots


## South Subcluster
#fitsfile = '/Users/dawson/Git/dissertation/Chapter4/BVRzMed_24maglim_pzpen_south_10000_numberdensity.fits'
#x_start = 148
#y_start = 50
#r_start = 40
#prefix = 'galden_pzpen_10000_south'

#cent_n_var(fitsfile,id_sci,N_boot,id_boot,x_start,y_start,r_start,r_step,r_stop,N_bins,prefix)

# North Subcluster
fitsfile = '/Users/dawson/Git/dissertation/Chapter4/BVRzMed_24maglim_pzpen_north_10000_numberdensity.fits'
x_start = 104
y_start = 163
r_start = 40
prefix = 'galden_pzpen_10000_north'

cent_n_var(fitsfile,id_sci,N_boot,id_boot,x_start,y_start,r_start,r_step,r_stop,N_bins,prefix)

pylab.show()
