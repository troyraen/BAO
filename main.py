import numpy as np
import pandas as pd
from MockBox import MockBox as MB
# import setup as su
# import calc_wtheta as cw
# import myplots as mp
# import helper_fncs as hf

# # NOTE:
# see http://gouthamanbalaraman.com/blog/numpy-vs-pandas-comparison.html
# for operation speed comparisons
#

cat_gals = 5e5 # approx num gals in cat mock
Nrands = 10*cat_gals # cat has ~5e5 galaxies
# tbin_edges = np.logspace(np.log10(0.1), np.log10(15.0), 75)
tbin_edges = np.linspace(0.1, 15.0, 75)
zbin_width = 0.4 # [0.3, 0.4, 0.5]
# for zw in zbin_width:
zw = zbin_width
imax=10
for i in range(imax):
    print('\n\n*** Starting interation {} out of {}\n\n'.format(i+1,imax))
    mb = MB(Nstack=2, zbin_width=zw, tbin_edges=tbin_edges, Nrands=Nrands, galplots=False)
    mb.getmock_calcwtheta()#fow='all')
