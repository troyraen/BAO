# %matplotlib qt
import numpy as np
import pandas as pd
from astropy import cosmology
from astropy.table import Table
import datetime
import time

import setup as su
import calc_wtheta as cw
import myplots as mp
import helper_fncs as hf


def getmock_calcwtheta(Nstack=2, z4push=su.catboxz, zspace=0.365, tbins=None, \
        fout='data/wtheta.dat', zrunfout='data/zruntime.dat', nthreads=32, galplots=False):
    """
    Stacks Nstack^3 boxes together (around the origin) to create a bigger box.
    Needs update so Nstack=0 => just move origin to center of box for push consistency.
    Pushes the box so the face is at comoving_distance(redshift = su.catboxz (loaded in load_popmock))
    Transforms to RA, DEC, Z and bins redshift using zspace.
    Calculates wtheta using tbins, nthreads and writes results to fout.
    Calculates runtime of each wtheta calculation and outputs info to zrunfout.
    """
    # Notes:
    # All globals are defined in setup.py
    # Expect BAO at ~6.6 degrees for z~0.5
    # Halotools assumes all lengths are in Mpc/h
    # zspace: max redshift error in SDSS DR10 Photoz table is 0.365106,
        # see http://skyserver.sdss.org/CasJobs/MyDB.aspx MyTable_1 and
        # http://skyserver.sdss.org/dr6/en/help/docs/algorithm.asp?key=photoz
    # zrunfout code to print header:
        # zrunhdr = ['nthreads', 'zbin', 'calctime', 'numgals']
        # zrunhdrstr = ''.join(str(x).rjust(16) for x in zrunhdr)
        # print(zrunhdrstr, file=open(zrunfout, 'a'))


    # Setup:
    start_script = time.time() # time the script
    print('\ndo_mock_wtheta.py started at {}'.format(datetime.datetime.now()))
    su.load_cosmo() # loads global cosmo object plus H0, Om0
    su.load_popmock() # unnecessary, called from get_galtbl
    galdf = su.get_galtbl(getas='DF') # get the galaxy_table as a DataFrame
    if galplots:
        mp.plot_galaxies(galdf, gal_frac=0.005, coords='xyz', title='Original Mock') # plot a random subsample of galaxies
    if tbins is None:
        tbins = np.logspace(np.log10(1.0), np.log10(10.0), 15)
    print('*** You should update do_mock_wtheta to check if zrunfout exists, create/write header if not. ***')
    print('\t\t*** getmock_calcwtheta has header code in Notes. ***')


    # Stack boxes, push to catalog redshift, and transform coordinates
    print('Stacking {}^3 boxes. ...'.format(Nstack))
    newgals = su.stack_boxes(galdf, Nstack=Nstack, ogLbox=su.catLbox) # returns df
    if galplots:
        # ngtbl = Table(newgals, names=['x','y','z','vx','vy','vz'])
        mp.plot_galaxies(newgals, gal_frac=5e-4, coords='xyz', title='Boxes Stacked Around Origin')

    print('Pushing the box out to box x-face redshift = {0:1.2f} ...'.format(z4push))
    newgals = su.push_box2z(newgals, z4push, su.newLbox) # returns original DF with 'x' column shifted to so box x-face is at redshift z4push
    if galplots:
        # ngtbl = Table(newgals_atz, names=['x','y','z','vx','vy','vz'])
        mp.plot_galaxies(newgals, gal_frac=5e-4, coords='xyz', title='Boxes Stacked and Pushed to Catalog Redshift')

    print('Converting to RA, DEC, Redshift. ...')
    newgals = hf.get_ra_dec_z(newgals, usevel=True) # Adds columns to newgals

    # Bin redshifts calculate wtheta for each zbin and write to file
    newgals, zbin_edges = hf.bin_redshifs(newgals, zspace=zspace, validate=False)
    print('*** You should fix redshift bins so you get consistent binning with different mocks. ***')
    print('\t\t*** do_mock_wtheta.py line 64. ***')
    # zbcens = rdz.zbin.unique() # get set of zbin centers to use as masks
    zgroups = newgals.groupby('zbin') # group by redshift bin
    randoms_kwargs = { 'boxsize':su.newLbox, 'push_to_z':su.catboxz }
    mocknum = get_mock_num() # get mock number as date and time
    for i, (zzz, rdz_z) in zgroups:
        print('\nCalculating wtheta for zbin = {0:1.2f}\n\t{1}\n'.format(zzz, datetime.datetime.now()))
        start_zbin = time.time() # time the wtheta calculation
        # rdz_z = rdz.loc[rdz.zbin == zzz]
        tbcens, wtheta = cw.calc_wtheta(rdz_z, tbins, randoms_kwargs, nthreads=nthreads)
        cw.write_to_file(tbcens, wtheta, zzz, mocknum, fout)

        # print/save calculation time
        end_zbin = time.time() # time the wtheta calculation
        ztime = (end_zbin-start_zbin)/60. # in minutes
        zrundat = np.asarray([nthreads, zzz, ztime, len(rdz_z.index)])
        zrunstr = np.array2string(zrundat, formatter={'float_kind':lambda x: "%15.1f" % x})[1:-1]
        print(zrunstr, file=open(zrunfout, 'a'))
        print('wtheta calculation took {0:.1f} minutes with nthreads = {1}\n'.format(ztime, nthreads))
        print('Results written to {0}. Calc time written to {1}.'.format(fout, zrunfout))



    end_script = time.time() # time the script
    stime = (end_script-start_script)/60. # in minutes
    print('\n\t{0}\ndo_mock_wtheta.py ran for {1:.1f} minutes.\n'.format(datetime.datetime.now(), stime))


def get_mock_num():
    dtm = datetime.datetime.now() # get date and time to use as mock number
    mocknum = float(dtm.strftime("%m%d%y.%H%M"))
    return mocknum
