%matplotlib qt
import numpy as np
from astropy import cosmology
from astropy.table import Table

import setup_mock as sm
import calc_wtheta as cw
import myplots as mp


halocat, HODmodel = sm.setup_halosHOD() # fetch halo catalog and HODmodel, returns populated HODmodel.mock
HODmodel.mock.populate() # repopulate
galaxy_table = HODmodel.mock.galaxy_table # get the galaxy_table
# mp.plot_galaxies(galaxy_table, gal_frac=0.005, coords='xyz') # plot a random subsample of galaxies

# stack Nstack boxes together to create a bigger box
boxsize = halocat.Lbox[0]
newgals = sm.stack_boxes(galaxy_table, Nstack=2, Lbox=boxsize)
ngtbl = Table(newgals, names=['x','y','z','vx','vy','vz'])
mp.plot_galaxies(ngtbl, gal_frac=5e-4, coords='xyz')

# push the box out to x -> x + comoving_distance(halocat redshift)
cosmo = cosmology.FlatLambdaCDM(H0=0.7, Om0=0.3)
zred = halocat.redshift
xz = (cosmo.comoving_distance(zred).value)*cosmo.h # Mpc/h
# ???? multiply or divide here? Duncan's ra_dec_z divides to get RID of h scaling?
newgals[:,0] = newgals[:,0]+ xz

rdz = sm.get_ra_dec_z(newgals)


# #### SAND ####
# gt = galaxy_table[1:5]
# newgals = sm.stack_boxes(gt, Nstack=4, Lbox=boxsize)
# mp.plot_galaxies(gt, gal_frac=1, coords='xyz')
# ng = Table(newgals, names=['x','y','z','vx','vy','vz'])
# #### SAND ####
#
# # get box as dataframe with ra, dec, z
# # slice the box in redshift bins
# # calculate and store wtheta for each zbin
#
# bins = np.logspace(np.log10(5.0), np.log10(15.0), 15)
# bcens, wtheta = cw.get_wtheta(halocat, HODmodel, bins, fout='wtheta.dat')
#
# mp.plot_wtheta(bcens, wtheta)
