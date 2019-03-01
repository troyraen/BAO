# %matplotlib qt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy import cosmology
from astropy.table import Table
import datetime
import time

# import setup_mock as sm
impoort do_mock_wtheta as dmw
import calc_wtheta as cw
import myplots as mp
# import helper_fncs as hf

# start_script = time.time() # time the script
# print('\nmain.py started at {}'.format(datetime.datetime.now()))

fout = 'wtheta.dat'
zrunfout='zruntime.dat'
nthreads = 32
# tbins = np.logspace(np.log10(4.0), np.log10(8.0), 20)
tbins = np.linspace(4.0, 8.0, 20)
# if you change tbins, move the existing file ### FIX THIS ### (bash code)
# !tm=$(date +"%m%d%y_%H%M")
# !mv data/wtheta.dat data/wtheta_ow_${tm}.dat
#
dmw.getmock_calcwtheta(Nstack=2, zspace=0.365, tbins=tbins, \
        fout=fout, zrunfout=zrunfout, nthreads=nthreads, galplots=False)

fin = fout
wdf = cw.load_from_file(fin)
# pd.pivot_table(wdf, index='zbin')

zrf = mp.getplot_zruntimes() # get and plot zrun calc times


# end_script = time.time() # time the script
# stime = (end_script-start_script)/60. # in minutes
# print('\n\t{0}\nmain.py ran for {0:.1f} minutes.\n'.format(datetime.datetime.now(), stime))
