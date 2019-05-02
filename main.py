import numpy as np
import pandas as pd
from MockBox import MockBox as MBox
# import setup as su
import calc_wtheta as cw
import myplots as mp
# import helper_fncs as hf

# # NOTE:
# see http://gouthamanbalaraman.com/blog/numpy-vs-pandas-comparison.html
# for operation speed comparisons
#

# DEFAULTS:
stats=['wtheta', 'xi', 'wp']
nbins = 50
tbin_edges = np.logspace(np.log10(0.1), np.log10(12.0), nbins+1)
rbin_edges = np.logspace(np.log10(100.0), np.log10(160.0), nbins+1)
pimax = 500
z4push = 'cat'
zw = 2.
cat_gals = 5e5 # approx num gals in cat mock
Nstack = 0 # number of mock boxes to stack, per dimension
nrfact = 3 # used in Nrands
Nrands = int(nrfact*cat_gals* max(Nstack,1)**3)
imax = 1 # number of times to run each param combo


for i in range(imax):
    print('\n\n*** Starting interation {} out of {}\n\n'.format(i+1,imax))
    mb = MBox(Nstack=Nstack, z4push=z4push, zbin_width=zw, tbin_edges=tbin_edges, \
                rbin_edges=rbin_edges, pimax=pimax, Nrands=Nrands, galplots=False)
    mb.getmock()
    mb.calc_stats(stats=stats, nthreads=24)




#
# Nstack_lst = [0]#, 2]
# nrfact_lst = [1]#,3,10]
# for Nstack in Nstack_lst:
#     print('\n-------------------------------------------------------------------')
#     print('\n*** Starting Nstack = {}\n\n'.format(Nstack))
#
#     for nrfact in nrfact_lst:
#         print('\n____________________________________________________________________')
#         print('\n*** Starting nrfact = {}\n\n'.format(nrfact))
#
#         Nrands = int(nrfact*cat_gals* max(Nstack,1)**3)
#         for i in range(imax):
#             print('\n\n*** Starting interation {} out of {}\n\n'.format(i+1,imax))
#             mb = MBox(Nstack=Nstack, z4push=z4push, zbin_width=zw, tbin_edges=tbin_edges, rbin_edges=rbin_edges, Nrands=Nrands, galplots=False)
#             # mb.getmock_calcwtheta(nthreads=24)
#             mb.getmock()
#             mb.calc_stats(stats=stats)


# # Get and plot wtheta from file
# fin = 'data/wtheta.dat'
# wdf = cw.load_from_file(fin)
# wdf_nstackg = wdf.groupby('Nstack')
# for i, (Ns, wdfn) in enumerate(wdf_nstackg):
#     mp.plot_wtheta(wdfn, save='plots/wtheta_Nstack{}.png'.format(Ns), show=False)
