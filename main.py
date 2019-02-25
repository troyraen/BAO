# %matplotlib qt
import numpy as np
import pandas as pd
from astropy import cosmology
from astropy.table import Table

import setup_mock as sm
import calc_wtheta as cw
import myplots as mp
import helper_fncs as hf


# Halotools assumes all lengths are in Mpc/h
H0 = 70.0
try:
    halocat and HODmodel
except:
    halocat, HODmodel = sm.setup_halosHOD() # fetch halo catalog and HODmodel, returns populated HODmodel.mock
    # HODmodel.mock.populate() # repopulate
    catLbox = halocat.Lbox[0]
    catboxz = halocat.redshift

galaxy_table = HODmodel.mock.galaxy_table # get the galaxy_table
# mp.plot_galaxies(galaxy_table, gal_frac=0.005, coords='xyz') # plot a random subsample of galaxies

# stack Nstack^3 boxes together to create a bigger box
Nstack = 2
newLbox = catLbox*Nstack
newgals = sm.stack_boxes(galaxy_table, Nstack=Nstack, Lbox=catLbox) # returns (ngals x 6) ndarray
# ngtbl = Table(newgals, names=['x','y','z','vx','vy','vz'])
# mp.plot_galaxies(ngtbl, gal_frac=5e-4, coords='xyz')

# push the box out to x -> x + comoving_distance(halocat redshift)
cosmo = cosmology.FlatLambdaCDM(H0=H0, Om0=0.3)
newgals_atz = sm.push_box2z(newgals, catboxz, newLbox, cosmo=cosmo) # returns original ndarray with 1st column shifted
# xzbox = (cosmo.comoving_distance(catboxz).value)*cosmo.h # Mpc/h
# newgals[:,0] = newgals[:,0]+ xzbox
# ngtbl = Table(newgals, names=['x','y','z','vx','vy','vz'])
# mp.plot_galaxies(ngtbl, gal_frac=5e-4, coords='xyz')

# transform to ra, dec, and redshift
# rdzF = pd.DataFrame(hf.get_ra_dec_z(newgals, cosmo=cosmo, usevel=False), columns=['RA','DEC','Redshift'])
# rdzT = pd.DataFrame(hf.get_ra_dec_z(newgals, cosmo=cosmo, usevel=True), columns=['RA','DEC','Redshift'])
# mp.plot_galaxies(rdz, gal_frac=5e-4, coords='rdz')
rdz = hf.get_ra_dec_z(newgals_atz, cosmo=cosmo, usevel=True) # now returns a df
# rdz = rdzT

# bin redshifts
# hf.bin_redshifs(rdz, zspace)
# for each mask, get slice of rdz and calc wtheta
# save results in df (and write to file) with
    # columns: {wtheta bcens (1 col each, holds wtheta for the bin),
    #           zbin (hold zbcens value for given row),
    #           mock number ()}
    # rows: {zbcen (1 row each)}

# then calculate wtheta
