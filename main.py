# %matplotlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy import cosmology
from astropy.table import Table
import datetime
import time

import setup as su
import do_mock_wtheta as dmw
import calc_wtheta as cw
import myplots as mp
# import helper_fncs as hf


# gen mock, calc wtheta, and write to file
fout = 'data/wtheta.dat'
zrunfout='data/zruntime.dat'
nthreads = 32
tbins = np.logspace(np.log10(1.0), np.log10(13.0), 50)
# tbins = np.linspace(4.0, 8.0, 20)
# z4push = 0.
Nstack=0
su.load_popmock()
for z4push in [0., su.catboxz]:
    dmw.getmock_calcwtheta(Nstack=Nstack, z4push=z4push, zspace=0.365, tbins=tbins, \
        fout=fout, zrunfout=zrunfout, nthreads=nthreads, galplots=False)
###
