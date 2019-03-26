###
### This is a class to handle the data associated with a mock box populated with galaxies.
###
### Troy Joseph Raen
### March 2019 (started)


import numpy as np
import pandas as pd
from collections import OrderedDict as OD
# from astropy import cosmology
# from astropy.table import Table
import datetime
# import time

import setup as su
import calc_wtheta as cw
import myplots as mp
import helper_fncs as hf



class MockBox:
    def __init__(self, Nstack=0, zbin_width=0.365, tbin_edges=None, rtfout='data/runtimes.dat'):
        self.mocknum = get_mock_num() # float. Date, time formatted as "%m%d%y.%H%M"
        self.report_times = None # ordered dict for fncs to report runtimes. Generally, key = fnc name, val = (start time while fnc running, overwrite with:) fnc runtime
        self.rtfout = rtfout # = string writes function runtimes to this file. = None skips timing fncs.

        self.Lbox = None # scalar. Length of box side (assumed square). Mpc/h
        self.zbox = None # scalar. Redshift at face of box.
        self.zbins = None # array. Redshift at center of each redshift bin.
        self.zbin_edges = None # array. Redshift bin edges.
        self.zbin_width = zbin_width # desired spacing between bin centers (result will not be exact)

        self.PhaseSpace = None # DataFrame. Columns {x,y,z, vx,vy,vz}, rows = galaxies.
        self.RDZ = None # DataFrame. Columns {RA, DEC, Redshift}, rows = galaxies. bin_redshifs() adds column zbin
        self.Randoms = None # DataFrame. Columns {RA, DEC, zbin}, rows = random points within a box defined by self.Lbox and self.zbox. zbin should be one of self.zbins
        self.tbins = None # array. Theta [deg] at center of bin.
        self.tbin_edges = tbin_edges # array. Theta bin edges. if == None, set default below
        self.wtheta = None # DataFrame. Columns {theta bins}, rows wtheta, index zbin

        # Data from original (halos + HOD galaxies) catalog mock
        self.cat_galtbl = None  # DataFrame. Columns {x,y,z, vx,vy,vz}, rows = galaxies. This should come from populating a mock with HOD.
        self.cat_Lbox = None # Length of mock box side (assumed square). Mpc/h
        self.cat_zbox = None # scalar. Redshift at center? of box.
        self.Nstack = Nstack # even integer. number of catalog mock boxes to stack together

        if tbin_edges is None:
            self.tbin_edges = np.logspace(np.log10(1.0), np.log10(10.0), 25)
            print('tbin_edges instantiated to: {}'.format(self.tbin_edges))


    # Methods:
    def get_mock_num(self):
        dtm = datetime.datetime.now() # get date and time to use as mock number
        mocknum = float(dtm.strftime("%m%d%y.%H%M"))
        return mocknum

    # push_box2catz = push_box2catz
    # stack_boxes = stack_boxes
    # stack_boxes_gen_new_rows = stack_boxes_gen_new_rows
    # transform_mock = transform_mock
    # find_bin_center = find_bin_center
    # bin_redshifs = bin_redshifs
    # calc_write_wtheta = calc_write_wtheta
    # getmock_calcwtheta = getmock_calcwtheta # generate new mock then calc wtheta



    def getmock_calcwtheta(self, tbin_edges=None, fout='data/wtheta.dat', Nstack=None, rtfout=None, zbin_width=None, nthreads=32, galplots=False):
        """
        Stacks Nstack^3 boxes together (around the origin) to create a bigger box.
        DONE: Needs update so Nstack=0 => just move origin to center of box for push consistency.
        Pushes the box so the face is at comoving_distance(redshift = self.cat_zbox)
        Transforms to RA, DEC, z and bins redshift using self.zbin_width.
        Calculates wtheta using tbins, nthreads and writes results to fout.
        Calculates runtime of each wtheta calculation and outputs info to rtfout.
            rtfout == string writes function runtimes to this file
                   == None to skip timing fncs.

        Notes:

        All globals are defined in setup.py
        Expect BAO at ~6.6 degrees for z~0.5
        Halotools assumes all lengths are in Mpc/h
        self.zbin_width: max redshift error in SDSS DR10 Photoz table is 0.365106,
            see http://skyserver.sdss.org/CasJobs/MyDB.aspx MyTable_1 and
            http://skyserver.sdss.org/dr6/en/help/docs/algorithm.asp?key=photoz
            Note from Jeff: SDSS has unusually large z errors.
        """


        ### Setup:
        if rtfout is not None: # this is set on __init__, but can be changed here
            self.rtfout = rtfout
        if self.rtfout is not None: # set up dict to track function runtimes
            self.report_times = OD([('mocknum', mocknum), ('nthreads',nthreads)]) # report_times dict for functions to report times
            self.report_times['fcol_width'] = 25 # set report_times file column width
            self.report_times['getmock_calcwtheta'] = hf.time_code('start') #.TS. get code start time
            print('Function runtimes will be written to {}'.format(self.rtfout))
            print('getmock_calcwtheta() started at {}'.format(datetime.datetime.now()))
        else:
            print('*** Function runtimes are NOT being tracked.\n\t Set MockBox.rtfout to track these.')

        if Nstack is not None: # this is set on __init__, but can be changed here
            self.Nstack = Nstack

        if zbin_width is not None: # this is set on __init__, but can be changed here
            self.zbin_width = zbin_width

        if tbin_edges is not None: # this is set on __init__, but can be changed here
            self.tbin_edges = tbin_edges

        su.load_cosmo() # loads default global cosmo object plus H0, Om0
        ###

        # Get galaxy DF by populating DM mock using HODmodel
        self.cat_galtbl, self.cat_Lbox, self.cat_zbox = su.get_galtbl(getas='DF')
        if galplots: # plot a random subsample of galaxies
            mp.plot_galaxies(self.cat_galtbl, gal_frac=0.005, coords='xyz', title='Original Mock')

        # Stack boxes and push box face to appropriate redshift.
        self.transform_mock()

        # Transform coordinates
        print('\nConverting to RA, DEC, Redshift. ...')
        if self.rtfout is not None:
            self.report_times['get_ra_dec_z'] = hf.time_code('start') #.TS. get code start time
        self.RDZ = hf.get_ra_dec_z(self.PhaseSpace, usevel=True)
        if self.rtfout is not None: # set up dict to track function runtimes
            self.report_times['get_ra_dec_z'] = hf.time_code(self.report_times['get_ra_dec_z'], unit='min') #.TE. replace start time with runtime in minutes


        # Bin redshifts, add column to self.RDZ
        self.bin_redshifs(validate=False)

        # Calculate wtheta for each zbin and write to file
        zgroups = self.RDZ[['RA','DEC','zbin']].groupby('zbin', axis=0) # group by redshift bin, only need these columns
        for i, (zzz, rdz_z) in enumerate(zgroups):
            # Several report_times entries are set in calc_write_wtheta and overwritten for each zzz.
            # Be sure to write report_times to file before the end of this for loop
            self.calc_write_wtheta(zzz, rdz_z, fout, nthreads=nthreads)

            if self.rtfout is not None:
                hf.write_report_times(self.report_times, rtfout)
                print('\tcalc_wtheta for zbin = {0} took {1:.1f} minutes'.format(zzz, self.report_times['calc_wtheta']))
                print('Results written to {0}. Calculation report_times written to {1}.\n'.format(fout, rtfout))
                # print('{0:15d} {1:15.1f} {2:15.1f} {3:15d} {4:15.4f} {5:15d}'.format(nthreads, zzz, ztime, len(rdz_z.index), dtm, len(tbins)), file=open(rtfout, 'a'))
                # print('\nwtheta calculation took {0:.1f} minutes with nthreads = {1}\n'.format(ztime, nthreads))


        if self.rtfout is not None:
            self.report_times['getmock_calcwtheta'] = hf.time_code(self.report_times['getmock_calcwtheta'], unit='min') #.TE. replace start time with runtime in minutes
            print('\n\t{0}\ndo_mock_wtheta.py ran for {1:.1f} minutes.\n'.format(datetime.datetime.now(), self.report_times['getmock_calcwtheta']))

        return None



    def calc_write_wtheta(self, zzz, rdz, fout, nthreads=32):

        if self.rtfout is not None: # Save some info to report_times
            self.report_times['zbin'] = zzz
            self.report_times['numgals'] = len(rdz.index)
            print('\nCalculating wtheta for zbin = {0:1.2f}\n\t{1}\n'.format(zzz, datetime.datetime.now()))
            self.report_times['calc_wtheta'] = hf.time_code('start') #.TS. get code start time

        randoms_kwargs = { 'Nran':len(rdz.index), 'viewgals':galplots }
        tbcens, wtheta, self.report_times = cw.calc_wtheta(rdz, MockBox=self, randoms_kwargs=randoms_kwargs, nthreads=nthreads)

        # Set self.tbins or check that it equals tbcens
        if self.tbins is None:
            self.tbins = tbcens
        else:
            errmsg = 'Theta bin centers from cw.calc_wtheta do not match those previously set in MockBox.tbins.'
            np.testing.assert_allclose(tbcens, self.tbins, rtol=1e-5, err_msg=errmsg)
            # stops execution with error if theta bin centers to not match

        # Set or append to self.wtheta
        tmpdf = pd.DataFrame(data=wtheta, index=zzz, columns=tbcens)
        self.wtheta = tmpdf if self.wtheta is None else self.wtheta.append(tmpdf, ignore_index=False, sort=True)

        # Write current zbin wtheta info to file
        # Do this now so it is not lost if a future zbin calculation fails
        cw.write_to_file(tbcens, wtheta, zzz, mocknum, fout)

        if self.rtfout is not None:
            self.report_times['calc_wtheta'] = hf.time_code(self.report_times['calc_wtheta'], unit='min') #.TE. replace start time with runtime in minutes

        return None



    # SDSS DR10 max photo zErr = 0.365106
    def bin_redshifs(self, zbin_width=None, validate=False):
        """ Adds column 'zbin' to self.RDZ containing center of zbin the galaxy is in
                (zbin is one of self.zbins)
            Sets self.zbins and self.zbin_edges

            zbin_width = desired spacing between bins (result will not be exact)
        """

        if self.rtfout is not None: # track function runtimes
                self.report_times['bin_redshifs'] = hf.time_code('start') #.TS. get code start time

        if zbin_width is not None: # this is set on __init__, but can be changed here
            self.zbin_width = zbin_width

        print('\n*** You should fix redshift bins so you get consistent binning with different mocks. ***')

        zmin = self.RDZ.Redshift.min()
        zmax = self.RDZ.Redshift.max()
        # eps = 0.01
        # self.zbin_edges = np.arange(zmin, zmax+eps, self.zbin_width)
        num_binedges = int(np.ceil((zmax-zmin)/self.zbin_width))
        self.zbin_edges = np.linspace(zmin,zmax,num_binedges)
        self.zbins = np.round((self.zbin_edges[:-1]+self.zbin_edges[1:])/2, 2) # keep 2 decimal places
        # create zbin masks for rdz dataframe
            # add a column to rdz containing the zbin (self.zbins value) the galaxy resides in
        # given z, find which bin it's in and get the value of the bin center
        self.RDZ['zbin'] = self.RDZ['Redshift'].apply(find_bin_center, **{"bin_edges":self.zbin_edges, "bin_centers":self.zbins})
        # rdz.apply(lambda inval: hf.find_bin_center(inval, self.zbin_edges, self.zbins), rdz.Redshift)

        # make sure the operation worked as expected:
        if validate:
            for index, row in self.RDZ.iterrows():
                rdzbin = row.zbin
                truezbin = find_bin_center(row.Redshift, bin_edges=self.zbin_edges, bin_centers=self.zbins)
                assert(rdzbin==truezbin)


        if self.rtfout is not None: # track function runtimes
                self.report_times['bin_redshifs'] = hf.time_code(self.report_times['bin_redshifs'], unit='min') #.TE. replace start time with runtime in minutes

        return None



    def find_bin_center(self, inval, bin_edges=None, bin_centers=None):
        """Returns the bin_centers value corresponding to the
            bin (defined by bin_edges) that inval (scalar) falls in to.
        """
        print('\n*** find_bin_center has not been updated since moving to MockBox class.\n')
        for i in range(len(bin_centers)):
            if (bin_edges[i] <= inval) and (inval < bin_edges[i+1]):
                return bin_centers[i]
        if inval == bin_edges[-1]: # inval == top edge of last bin
            return bin_centers[-1]
        # if you get here, inval is not in any of the bins
        print('{} did not fall within any bins.'.format(inval))
        return None



    def transform_mock(self):
        """ Stack boxes around origin.
            Push face to z.
            Set self.PhaseSpace.

            All self.mock_* should be set prior to calling this.
                e.g. self.cat_galtbl, self.cat_Lbox, self.cat_zbox = su.get_galtbl(getas='DF')
         """

        # Stack boxes
        self.stack_boxes() # Sets self.PhaseSpace by stacking cat_galtbl, centered on origin.
        if galplots:
            mp.plot_galaxies(self.PhaseSpace, gal_frac=5e-4, coords='xyz', title='Boxes Stacked Around Origin')

        # Push box face to redshift = self.cat_zbox
        self.push_box2catz()
        if galplots:
            mp.plot_galaxies(self.PhaseSpace, gal_frac=5e-4, coords='xyz', title='Boxes Stacked and Pushed to Redshift = {}'.format(self.zbox))


    def stack_boxes_gen_new_rows(self, gal):#, Nstack=None, ogLbox=None):
        """ This gets applied to each row of the original dataframe.
            gal = single galaxy data as series (single row from original df)
            Generates new 'xyz' coordinates by stacking self.Nstack**3 boxes.
            Copy all other columns.
            Returns a DataFrame with self.Nstack**3 galaxies (rows)
        """
        N2 = int(self.Nstack/2)
        assert Nstack%2 == 0, 'Nstack should be an even integer.'

        dfcols = list(gal.index.values)
        extra_vals = gal[[c for c in dfcols if c not in ['x','y','z']]] # get series of just extra values

        numcols = len(dfcols)
        gtot = self.Nstack**3 # num galaxies generated from single original galaxy
        newgals = np.zeros((gtot,numcols)) # create an array to hold the output
        nLin = np.linspace(-N2,N2-1,self.Nstack)*self.cat_Lbox  # array of ints -N2 to N2-1
        # boxes will be stacked around the origin by adding nLin to each galaxy coordinate
        # create coordinates for stacked boxes in each direction
        xstk = gal.x*np.ones(self.Nstack) +nLin # array of length Nstack
        ystk = gal.y*np.ones(self.Nstack) +nLin
        zstk = gal.z*np.ones(self.Nstack) +nLin
        gnew = 0
        for i in xstk:
            for j in ystk:
                for k in zstk:
                    coords = np.array([ i, j, k ] + list(extra_vals))
                    newgals[gnew] = coords # add the new galaxy
                    gnew = gnew+1
        ngdf = pd.DataFrame(newgals, columns=dfcols)
        return ngdf



    def stack_boxes(self):#galdf, Nstack=2, ogLbox=1000):
        """
        galdf = DataFrame representing a mock box,
                with minimum columns {'x','y','z', 'vx','vy','vz'},
                (assumes coordinates {x,y,z} in [0,Lbox]).
        Stacks Nstack^3 (must be an EVEN integer) galdf mock boxes together around the origin.
        Sets self.Lbox.
        Returns a DataFrame of galaxies with the origin at the center of the box and boxsize=self.Lbox.
                columns {x,y,z} generated periodically, {vx,vy,vz, all others} same as original galaxies.
        """
        # Setup
        if self.rtfout is not None: # track function runtimes
            self.report_times['stack_boxes'] = hf.time_code('start') #.TS. get code start time

        self.Lbox = self.cat_Lbox*self.Nstack
        print('\nStacking {}^3 boxes. ...'.format(self.Nstack))
        print('*** Warning: MockBox.stack_boxes assumes the original box is strictly in the 1st quadrant with the origin at the corner. ***')

        if Nstack == 0: # Just move the origin to the center of the box.
            print('Moving origin to box center...')
            self.Lbox = self.cat_Lbox
            L2 = self.Lbox/2.
            newgals = self.cat_galtbl.copy(deep=True)
            newgals.x = newgals.x-L2; newgals.y = newgals.y-L2; newgals.z = newgals.z-L2

        else:
            # Iterate over rows of self.cat_galtbl.
            # Generate new 'xyz' coordinates. Copy all other columns.
            nglist = [] # list to hold each new DataFrame generated from a single galaxy
            for idx, gal in self.cat_galtbl.iterrows():
                nglist.append(stack_boxes_gen_new_rows(gal))
            # Create new df from dfs in nglist
            newgals = pd.concat(nglist, ignore_index=True)

        self.PhaseSpace = newgals

        if self.rtfout is not None: # track function runtimes
            self.report_times['stack_boxes'] = hf.time_code(self.report_times['stack_boxes'], unit='min') #.TE. replace start time with runtime in minutes
            print('\n\t{0}\nstack_boxes() ran for {1:.1f} minutes.\n'.format(datetime.datetime.now(), self.report_times['stack_boxes']))

        return None



    # cosmo.comoving_distance(zred) returns "Comoving distance in Mpc to each input redshift.""
    # from help(cosmo): Dimensionless Hubble constant: h = H_0 / 100 [km/sec/Mpc]
    def push_box2catz(self):#galdf, redshift, Lbox):
        """ self.PhaseSpace must be set before calling this function.
            Moves the self.PhaseSpace coordinates so the x-FACE is at
                comoving_distance(self.cat_zbox). (Only 'x' column is changed.)
            Sets self.zbox = self.cat_zbox
        """
        global cosmo
        if cosmo is None:
            load_cosmo()

        print('\nPushing the box out to box x-face redshift = {0:1.2f} ...'.format(self.cat_zbox))
        if self.rtfout is not None: # track function runtimes
            self.report_times['push_box2catz'] = hf.time_code('start') #.TS. get code start time

        # Shift x so coordinates are strictly positive (i.e. move face to x=0)
        # Then push the face to comoving_distance(self.cat_zbox)
        deltax = self.Lbox/2. + (cosmo.comoving_distance(self.cat_zbox).value)*cosmo.h # Mpc/h
        self.PhaseSpace.x = self.PhaseSpace.x + deltax

        self.zbox = self.cat_zbox

        if self.rtfout is not None: # track function runtimes
            self.report_times['push_box2catz'] = hf.time_code(self.report_times['push_box2catz'], unit='min') #.TE. replace start time with runtime in minute

        return None




# def NOTES_notrealfnc():
#
#     def myfunc(self):
#         print("Hello my name is " + self.name)
#
#     # multi indexing:
#     https://pandas.pydata.org/pandas-docs/stable/user_guide/advanced.html
#
#     # define method outside class:
#     # option 1:
#     def func(self):
#         print "func"
#
#     class MyClass(object):
#         myMethod = func
#
#     # option 2: after declaration
#     class MyClass(object):
#         pass
#
#     def func(self):
#         print "func"
#
#     MyClass.myMethod = func
