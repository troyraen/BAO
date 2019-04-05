# %matplotlib
import numpy as np
import pandas as pd
# import matplotlib.pyplot as plt
# import matplotlib as mpl
# from astropy import cosmology
# from astropy.table import Table
# import datetime
# import time

from MockBox import MockBox as MB
# import setup as su
# import do_mock_wtheta as dmw
# import calc_wtheta as cw
# import myplots as mp
# import helper_fncs as hf

# # NOTE:
# see http://gouthamanbalaraman.com/blog/numpy-vs-pandas-comparison.html
# for operation speed comparisons
#

mb = MB()
mb.getmock_calcwtheta(galplots=True, Nstack=2)




### OLD ###

# # gen mock, calc wtheta, and write to file
# fout = 'data/wtheta.dat'
# zrunfout='data/zruntime.dat'
# nthreads = 32
# tbins = np.logspace(np.log10(1.0), np.log10(13.0), 50)
# # tbins = np.linspace(4.0, 8.0, 20)
# su.load_popmock()
# z4push = su.catboxz
# galplots=False
# zspace=0.365
# Nstack=2
# # for z4push in [0., su.catboxz]:
# dmw.getmock_calcwtheta(Nstack=Nstack, z4push=z4push, zspace=zspace, tbins=tbins, \
#         fout=fout, zrunfout=zrunfout, nthreads=nthreads, galplots=galplots)
# ###
