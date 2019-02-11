import setup_halos as sh
import calc_wtheta as cw
import numpy as np
# %matplotlib qt
# import myplots as mp

halocat, HODmodel = sh.setup_halos()
bins = np.logspace(np.log10(5.0), np.log10(15.0), 15)
bcens, wtheta = cw.get_wtheta(halocat, HODmodel, bins)
