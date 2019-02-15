%matplotlib qt
import numpy as np

import setup_mock as sm
import calc_wtheta as cw
import myplots as mp


halocat, HODmodel = sm.setup_halosHOD() # fetch halo catalog and HODmodel, returns populated HODmodel.mock
HODmodel.mock.populate() # repopulate
galaxy_table = HODmodel.mock.galaxy_table # get the galaxy_table
mp.plot_galaxies(galaxy_table, gal_frac=0.005, coords='radecz') # plot a random subsample of galaxies
sm.stack_boxes(galaxy_table, Nstack=20)
# stack N sim boxes together to create a bigger box
# get box as dataframe with ra, dec, z
# slice the box in redshift bins
# calculate and store wtheta for each zbin

bins = np.logspace(np.log10(5.0), np.log10(15.0), 15)
bcens, wtheta = cw.get_wtheta(halocat, HODmodel, bins, fout='wtheta.dat')

mp.plot_wtheta(bcens, wtheta)
