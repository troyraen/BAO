# %matplotlib qt
import numpy as np
import pandas as pd
from astropy import cosmology
from astropy.table import Table
import datetime
import time

import setup_mock as sm
import calc_wtheta as cw
import myplots as mp
import helper_fncs as hf

start_script = time.time() # time the script
print('\nmain.py started at {}'.format(datetime.datetime.now()))

# Halotools assumes all lengths are in Mpc/h
H0 = 70.0
try:
    halocat and HODmodel
except:
    print('\nSetting up halocat and HODmodel\n')
    halocat, HODmodel = sm.setup_halosHOD() # fetch halo catalog and HODmodel, returns populated HODmodel.mock
    # HODmodel.mock.populate() # repopulate
    catLbox = halocat.Lbox[0]
    catboxz = halocat.redshift

galaxy_table = HODmodel.mock.galaxy_table # get the galaxy_table
# mp.plot_galaxies(galaxy_table, gal_frac=0.005, coords='xyz') # plot a random subsample of galaxies

# stack Nstack^3 boxes together to create a bigger box
Nstack = 2
newLbox = catLbox*Nstack
print('\nStacking {}^3 boxes\n'.format(Nstack))
newgals = sm.stack_boxes(galaxy_table, Nstack=Nstack, Lbox=catLbox) # returns (ngals x 6) ndarray
# ngtbl = Table(newgals, names=['x','y','z','vx','vy','vz'])
# mp.plot_galaxies(ngtbl, gal_frac=5e-4, coords='xyz')

# push the box out to x -> x + comoving_distance(halocat redshift)
cosmo = cosmology.FlatLambdaCDM(H0=H0, Om0=0.3)
print('\nPushing the box out to z={}\n'.format(catboxz))
newgals_atz = sm.push_box2z(newgals, catboxz, newLbox, cosmo=cosmo) # returns original ndarray with 1st column shifted
# xzbox = (cosmo.comoving_distance(catboxz).value)*cosmo.h # Mpc/h
# newgals[:,0] = newgals[:,0]+ xzbox
ngtbl = Table(newgals_atz, names=['x','y','z','vx','vy','vz'])
mp.plot_galaxies(ngtbl, gal_frac=5e-4, coords='xyz', title='Mock Galaxies')

# transform to ra, dec, and redshift
# rdzF = pd.DataFrame(hf.get_ra_dec_z(newgals, cosmo=cosmo, usevel=False), columns=['RA','DEC','Redshift'])
# rdzT = pd.DataFrame(hf.get_ra_dec_z(newgals, cosmo=cosmo, usevel=True), columns=['RA','DEC','Redshift'])
# mp.plot_galaxies(rdz, gal_frac=5e-4, coords='rdz')
print('\nConverting to RA, DEC, z\n')
rdz = hf.get_ra_dec_z(newgals_atz, cosmo=cosmo, usevel=True) # now returns a df
# rdz = rdzT

# bin redshifts
zspace = 0.365 # max redshift error in SDSS DR10 Photoz table is 0.365106
rdz, zbin_edges = hf.bin_redshifs(rdz, zspace=zspace, validate=False)
print('\nYou should fix redshift bins so you get consistent binning with different mocks.\n')
# get set of zbin centers to use as masks
zbcens = rdz.zbin.unique()
# calculate wtheta for each zbin (expect BAO at ~6.6 degrees for z~0.5)
tbins = np.logspace(np.log10(1.0), np.log10(10.0), 15)
randoms_kwargs = { 'boxsize':newLbox, 'push_to_z':catboxz, 'cosmo':cosmo }
fout = 'wtheta.dat'
dtm = datetime.datetime.now() # get date and time to use as mock number
mocknum = float(dtm.strftime("%m%d%y.%H%M"))
nthreads = 24
for zzz in zbcens:
    print('\nCalculating wtheta for zbin = {0}\n\t{1}\n'.format(zzz, datetime.datetime.now()))
    start_zbin = time.time() # time the wtheta calculation
    rdz_z = rdz.loc[rdz.zbin == zzz]
    tbcens, wtheta = cw.calc_wtheta(rdz_z, tbins, randoms_kwargs, nthreads=nthreads)
    cw.write_to_file(tbcens, wtheta, zzz, mocknum, fout)
    end_zbin = time.time() # time the wtheta calculation
    print('wtheta calculation took {0:.4f} seconds with nthreads = {1}\n'.format(end_zbin-start_zbin, nthreads))

wdf = cw.load_from_file(fout)

end_script = time.time() # time the script
print('\n\t{0}\nmain.py ran for {1:.4f} seconds.\n'.format(datetime.datetime.now(), end_script-start_script))



# for each mask, get slice of rdz and calc wtheta
# save results in df (and write to file) with
    # columns: {wtheta bcens (1 col each, holds wtheta for the bin),
    #           zbin (hold zbcens value for given row),
    #           mock number ()}
    # rows: {zbcen (1 row each)}

# then calculate wtheta
