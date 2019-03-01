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
import do_mock_wtheta as dmw
import calc_wtheta as cw
import myplots as mp
# import helper_fncs as hf


# gen mock, calc wtheta, and write to file
fout = 'data/wtheta.dat'
zrunfout='data/zruntime.dat'
nthreads = 32
# tbins = np.logspace(np.log10(4.0), np.log10(8.0), 20)
tbins = np.linspace(4.0, 8.0, 20) # if you change tbins, move the existing file ### FIX THIS ### (bash code)
# !tm=$(date +"%m%d%y_%H%M")
# !mv data/wtheta.dat data/wtheta_ow_${tm}.dat
#
dmw.getmock_calcwtheta(Nstack=2, zspace=0.365, tbins=tbins, \
        fout=fout, zrunfout=zrunfout, nthreads=nthreads, galplots=False)
###


# get wtheta df from file
import pandas as pd
import calc_wtheta as cw
fin = 'wtheta.dat'
wdf = cw.load_from_file(fin)
wdfp = pd.pivot_table(wdf, index='zbin')
###

# plot zbin wtheta calc times
import myplots as mp
zrf = mp.getplot_zruntimes(zrunfout='data/zruntime.dat') # get and plot zrun calc times
###
