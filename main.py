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
tbins = np.logspace(np.log10(3.0), np.log10(15.0), 20)
# tbins = np.linspace(4.0, 8.0, 20)
z4push = 0.
dmw.getmock_calcwtheta(Nstack=2, z4push=z4push, zspace=0.365, tbins=tbins, \
        fout=fout, zrunfout=zrunfout, nthreads=nthreads, galplots=True)
###
