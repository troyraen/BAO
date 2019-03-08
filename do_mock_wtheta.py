# %matplotlib qt
import numpy as np
import pandas as pd
from collections import OrderedDict as OD
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
    Pushes the box so the face is at comoving_distance(redshift = z4push)
    Transforms to RA, DEC, Z and bins redshift using zspace.
    Calculates wtheta using tbins, nthreads and writes results to fout.
    Calculates runtime of each wtheta calculation and outputs info to zrunfout.

    Notes:

    All globals are defined in setup.py
    Expect BAO at ~6.6 degrees for z~0.5
    Halotools assumes all lengths are in Mpc/h
    zspace: max redshift error in SDSS DR10 Photoz table is 0.365106,
        see http://skyserver.sdss.org/CasJobs/MyDB.aspx MyTable_1 and
        http://skyserver.sdss.org/dr6/en/help/docs/algorithm.asp?key=photoz
    zrunfout code to print header:
zrunhdr = ['nthreads', 'zbin', 'calctime', 'numgals', 'datetime', 'numtbins']
zrunhdrstr = ''.join(str(x).rjust(16) for x in zrunhdr)
print(zrunhdrstr, file=open(zrunfout, 'a'))
    """


    # Setup:
    mocknum = get_mock_num() # get mock number as date and time
    # rt['mocknum'] = mocknum # keep track of this
    rt = OD([('mocknum', mocknum), ('nthreads',nthreads)]) # report_times dict for functions to report times
    rt['fcol_width'] = 25 # set report_times file column width
    rt['getmock_calcwtheta'] = hf.time_code('start') #.TS. get code start time
    print('\ngetmock_calcwtheta() started at {}'.format(datetime.datetime.now()))

    su.load_cosmo() # loads global cosmo object plus H0, Om0
    su.load_popmock() # unnecessary, called from get_galtbl
    galdf = su.get_galtbl(getas='DF') # get the galaxy_table as a DataFrame
    if galplots:
        mp.plot_galaxies(galdf, gal_frac=0.005, coords='xyz', title='Original Mock') # plot a random subsample of galaxies
    if tbins is None:
        tbins = np.logspace(np.log10(1.0), np.log10(10.0), 25)
    print('*** You should update do_mock_wtheta to check if zrunfout exists, create/write header if not. ***')
    print('\t\t*** getmock_calcwtheta has header code in Notes. ***')
    ###

    # Stack boxes
    print('\nStacking {}^3 boxes. ...'.format(Nstack))
    rt['stack_boxes'] = hf.time_code('start') #.TS. get code start time
    newgals = su.stack_boxes(galdf, Nstack=Nstack, ogLbox=su.catLbox)  #.TC. Returns new DF with boxes stacked around the origin.
    rt['stack_boxes'] = hf.time_code(rt['stack_boxes'], unit='min') #.TE. replace start time with runtime in minutes
    print('\n\t{0}\nstack_boxes() ran for {1:.1f} minutes.\n'.format(datetime.datetime.now(), rt['stack_boxes']))
    if galplots:
        mp.plot_galaxies(newgals, gal_frac=5e-4, coords='xyz', title='Boxes Stacked Around Origin')

    # Push to redshift z4push
    print('\nPushing the box out to box x-face redshift = {0:1.2f} ...'.format(z4push))
    rt['push_box2z'] = hf.time_code('start') #.TS. get code start time
    newgals = su.push_box2z(newgals, z4push, su.newLbox) #.TC. returns original DF with 'x' column shifted to so box x-face is at redshift z4push
    rt['push_box2z'] = hf.time_code(rt['push_box2z'], unit='min') #.TE. replace start time with runtime in minute
    if galplots:
        mp.plot_galaxies(newgals, gal_frac=5e-4, coords='xyz', title='Boxes Stacked and Pushed to Redshift = {}'.format(z4push))

    # Transform coordinates
    print('\nConverting to RA, DEC, Redshift. ...')
    rt['get_ra_dec_z'] = hf.time_code('start') #.TS. get code start time
    newgals = hf.get_ra_dec_z(newgals, usevel=True) #.TC. Adds columns to newgals
    rt['get_ra_dec_z'] = hf.time_code(rt['get_ra_dec_z'], unit='min') #.TE. replace start time with runtime in minutes


    # Bin redshifts calculate wtheta for each zbin and write to file
    rt['bin_redshifs'] = hf.time_code('start') #.TS. get code start time
    newgals, zbin_edges = hf.bin_redshifs(newgals, zspace=zspace, validate=False) #.TC. "
    rt['bin_redshifs'] = hf.time_code(rt['bin_redshifs'], unit='min') #.TE. replace start time with runtime in minutes
    print('\n*** You should fix redshift bins so you get consistent binning with different mocks. ***')
    print('\t\t*** do_mock_wtheta.py line 64. ***')
    zgroups = newgals[['RA','DEC','zbin']].groupby('zbin', axis=0) # group by redshift bin, only need these columns
    randoms_kwargs = { 'boxsize':su.newLbox, 'push_to_z':z4push, 'viewgals':galplots }
    for i, (zzz, rdz_z) in enumerate(zgroups):
        # Setup. rt entries get set here and overwritten for each zzz.
        # Be sure to write to file before the end of this for loop
        rt['zbin'] = zzz # save this
        rt['numgals'] = len(rdz_z.index) # save this
        print('\nCalculating wtheta for zbin = {0:1.2f}\n\t{1}\n'.format(zzz, datetime.datetime.now()))
        rt['calc_wtheta'] = hf.time_code('start') #.TS. get code start time
        tbcens, wtheta, rt = cw.calc_wtheta(rdz_z, tbins, randoms_kwargs, nthreads=nthreads, report_times=rt) #.TC.
        rt['calc_wtheta'] = hf.time_code(rt['calc_wtheta'], unit='min') #.TE. replace start time with runtime in minutes
        cw.write_to_file(tbcens, wtheta, zzz, mocknum, fout)

        hf.write_report_times(rt, zrunfout)
        print('\tcalc_wtheta for zbin = {0} took {1:.1f} minutes'.format(zzz, rt['calc_wtheta']))
        print('Results written to {0}. Calculation report_times written to {1}.\n'.format(fout, zrunfout))
        # print('{0:15d} {1:15.1f} {2:15.1f} {3:15d} {4:15.4f} {5:15d}'.format(nthreads, zzz, ztime, len(rdz_z.index), dtm, len(tbins)), file=open(zrunfout, 'a'))
        # print('\nwtheta calculation took {0:.1f} minutes with nthreads = {1}\n'.format(ztime, nthreads))



    rt['getmock_calcwtheta'] = hf.time_code(rt['getmock_calcwtheta'], unit='min') #.TE. replace start time with runtime in minutes
    print('\n\t{0}\ndo_mock_wtheta.py ran for {1:.1f} minutes.\n'.format(datetime.datetime.now(), rt['getmock_calcwtheta']))

    return


def get_mock_num():
    dtm = datetime.datetime.now() # get date and time to use as mock number
    mocknum = float(dtm.strftime("%m%d%y.%H%M"))
    return mocknum
