import numpy as np
import pandas as pd
from MockBox import MockBox as MB
# import setup as su
import calc_wtheta as cw
import myplots as mp
# import helper_fncs as hf

# # NOTE:
# see http://gouthamanbalaraman.com/blog/numpy-vs-pandas-comparison.html
# for operation speed comparisons
#


tbin_edges = np.logspace(np.log10(0.1), np.log10(11.0), 76)
# tbin_edges = np.linspace(0.1, 6.0, 51)
z4push = 0.0
zw = 1.0 # 0.4 # [0.3, 0.4, 0.5] # zbin width
cat_gals = 5e5 # approx num gals in cat mock
Nstack_lst = [0]#, 2]
nrfact_lst = [1]#,3,10]

for Nstack in Nstack_lst:
    print('\n-------------------------------------------------------------------')
    print('\n*** Starting Nstack = {}\n\n'.format(Nstack))

    for nrfact in nrfact_lst:
        print('\n____________________________________________________________________')
        print('\n*** Starting nrfact = {}\n\n'.format(nrfact))

        Nrands = int(nrfact*cat_gals* max(Nstack,1)**3)
        imax=10
        for i in range(imax):
            print('\n\n*** Starting interation {} out of {}\n\n'.format(i+1,imax))
            mb = MB(Nstack=Nstack, z4push=z4push, zbin_width=zw, tbin_edges=tbin_edges, Nrands=Nrands, galplots=False)
            mb.getmock_calcwtheta(nthreads=24)


# get and plot wtheta from file
# fin = 'data/wtheta.dat'
# wdf = cw.load_from_file(fin)
# wdf_nstackg = wdf.groupby('Nstack')
# for i, (Ns, wdfn) in enumerate(wdf_nstackg):
#     mp.plot_wtheta(wdfn, save='plots/wtheta_Nstack{}.png'.format(Ns), show=False)
