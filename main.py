import numpy as np
import pandas as pd
from MockBox import MockBox as MBox
# import setup as su
# import calc_wtheta as cw
import myplots as mp
# import helper_fncs as hf



# DEFAULTS:
statfout='data/stats.dat'
stats=['wtheta', 'xi', 'wp']
nbins = 75
# tbin_edges = np.logspace(np.log10(0.05), np.log10(9.0), nbins+1)
tratio_binedges = np.logspace(np.log10(0.5), np.log10(1.5), nbins+1)
rbin_edges = np.logspace(np.log10(25.0), np.log10(150.0), nbins+1)
pimax = 300
galplots = False
z4push = 'cat'
zw_list = [0.1, 0.05, 0.15]
cat_gals = 5e5 # approx num gals in cat mock
Nstack_list = [0, 2] # number of mock boxes to stack, per dimension
nrfact_list = [1, 3] # used in Nrands
# Nrands = int(nrfact*cat_gals* max(Nstack,1)**3)
imax = 5 # number of times to run each param combo

for zw in zw_list:
    statfout='data/stats_tratiobins_zw{}.dat'.format(zw)
    for Nstack in Nstack_list:
        for nrfact in nrfact_list:

            Nrands = int(nrfact*cat_gals* max(Nstack,1)**3)
            print('----------------------------------------------------------------')
            print('Starting statfout = {}, Nstack = {}, nrfact = {}'.format(\
                            statfout, Nstack, nrfact))
            print('----------------------------------------------------------------')

            for i in range(imax):
                print('\n\n*** Starting interation {} out of {}\n\n'.format(i+1,imax))
                mb = MBox( Nstack=Nstack, z4push=z4push, zbin_width=zw, \
                            tratio_binedges=tratio_binedges, rbin_edges=rbin_edges, \
                            pimax=pimax, Nrands=Nrands, galplots=galplots, \
                            statfout=statfout
                        )
                mb.getmock()
                mb.calc_stats(stats=stats, nthreads=24)



## Do runs, Done
# fs

# # DEFAULTS:
# statfout='data/stats.dat'
# stats=['wtheta', 'xi', 'wp']
# nbins = 75
# tbin_edges = np.logspace(np.log10(0.05), np.log10(12.0), nbins+1)
# rbin_edges = np.logspace(np.log10(25.0), np.log10(150.0), nbins+1)
# pimax = 300
# galplots = False
# z4push = 'cat'
# zw = 10.
# cat_gals = 5e5 # approx num gals in cat mock
# Nstack = 2 # number of mock boxes to stack, per dimension
# nrfact = 1 # used in Nrands
# Nrands = int(nrfact*cat_gals* max(Nstack,1)**3)
# imax = 10 # number of times to run each param combo
#
# statfout='data/stats_Nstack0_z4push0.dat'
# print('----------------------------------------------------------------')
# print('Starting statfout = {}'.format(statfout))
# print('----------------------------------------------------------------')
# for i in range(imax):
#     print('\n\n*** Starting interation {} out of {}\n\n'.format(i+1,imax))
#     mb = MBox( Nstack=0, z4push=0., zbin_width=zw, tbin_edges=tbin_edges, \
#                 rbin_edges=rbin_edges, pimax=pimax, Nrands=Nrands, galplots=galplots, \
#                 statfout=statfout
#             )
#     mb.getmock()
#     mb.calc_stats(stats=stats, nthreads=24)
#
#
# statfout='data/stats_Nstack0.dat'
# print('----------------------------------------------------------------')
# print('Starting statfout = {}'.format(statfout))
# print('----------------------------------------------------------------')
# for i in range(imax):
#     print('\n\n*** Starting interation {} out of {}\n\n'.format(i+1,imax))
#     mb = MBox( Nstack=0, z4push=z4push, zbin_width=zw, tbin_edges=tbin_edges, \
#                 rbin_edges=rbin_edges, pimax=pimax, Nrands=Nrands, galplots=galplots, \
#                 statfout=statfout
#             )
#     mb.getmock()
#     mb.calc_stats(stats=stats, nthreads=24)
#
#
# statfout='data/stats_defaults.dat'
# print('----------------------------------------------------------------')
# print('Starting statfout = {}'.format(statfout))
# print('----------------------------------------------------------------')
# for i in range(imax):
#     print('\n\n*** Starting interation {} out of {}\n\n'.format(i+1,imax))
#     mb = MBox( Nstack=Nstack, z4push=z4push, zbin_width=zw, tbin_edges=tbin_edges, \
#                 rbin_edges=rbin_edges, pimax=pimax, Nrands=Nrands, galplots=galplots, \
#                 statfout=statfout
#             )
#     mb.getmock()
#     mb.calc_stats(stats=stats, nthreads=24)

### Older:
# statfout='data/stats_z4push0.dat'
# print('----------------------------------------------------------------')
# print('Starting statfout = {}'.format(statfout))
# print('----------------------------------------------------------------')
# for i in range(imax):
#     print('\n\n*** Starting interation {} out of {}\n\n'.format(i+1,imax))
#     mb = MBox( Nstack=Nstack, z4push=0., zbin_width=zw, tbin_edges=tbin_edges, \
#                 rbin_edges=rbin_edges, pimax=pimax, Nrands=Nrands, galplots=galplots, \
#                 statfout=statfout
#             )
#     mb.getmock()
#     mb.calc_stats(stats=stats, nthreads=24)
#



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
#             mb = MBox(Nstack=Nstack, z4push=z4push, zbin_width=zw, tbin_edges=tbin_edges, \
                    # rbin_edges=rbin_edges, Nrands=Nrands, galplots=False)
#             # mb.getmock_calcwtheta(nthreads=24)
#             mb.getmock()
#             mb.calc_stats(stats=stats)

# fe ## Do runs, Done

# # Get and plot wtheta from file
# fin = 'data/wtheta.dat'
# wdf = cw.load_from_file(fin)
# wdf_nstackg = wdf.groupby('Nstack')
# for i, (Ns, wdfn) in enumerate(wdf_nstackg):
#     mp.plot_wtheta(wdfn, save='plots/wtheta_Nstack{}.png'.format(Ns), show=False)
