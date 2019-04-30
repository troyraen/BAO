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
from pathlib import Path

import setup as su
import calc_wtheta as cw
import myplots as mp
import helper_fncs as hf



class MockBox:
    def __init__(self, Nstack=0, zbin_width=0.365, tbin_edges=None, wtfout='data/wtheta.dat', rtfout='data/runtimes.dat', Nrands=None, z4push='cat', galplots=False):
        self.mocknum = self.get_mock_num() # float. Date, time formatted as "%m%d%y.%H%M"
        self.report_times = OD([]) # ordered dict for fncs to report runtimes. Generally, key = fnc name, val = (start time while fnc running, overwrite with:) fnc runtime
        self.rtfout = rtfout # = string writes function runtimes to this file. = None skips timing fncs.
        self.galplots = galplots # bool. Whether to plot galaxy positions at each transformation (use to check whether transformations are correct.)

        # Data from original (halos + HOD galaxies) catalog mock
        self.cat_galtbl = None  # DataFrame. Columns {x,y,z, vx,vy,vz}, rows = galaxies. This should come from populating a mock with HOD.
        self.cat_Lbox = None # Length of mock box side (assumed square). Mpc/h
        self.cat_zbox = None # scalar. Redshift at center? of box.
        self.Nstack = Nstack # even integer. number of catalog mock boxes to stack together

        # Mock Box (stacked, coordinate transformed using original catalog mock)
        self.Lbox = None # scalar. Length of box side (assumed square). Mpc/h
        self.zbox = None # scalar. Redshift at face of box.
        self.z4push = z4push # = float. Redshift to push the faces of the boxes to (RDZ and Randoms).
                             # = 'cat' sets this to self.cat_zbox.
        self.zbins = None # array. Redshift at center of each redshift bin.
        self.zbin_edges = None # array. Redshift bin edges.
        self.zbin_width = zbin_width # desired spacing between bin centers (result will not be exact)
        self.PhaseSpace = None # DataFrame. Columns {x,y,z, vx,vy,vz, zbin}, rows = galaxies.
        self.RDZ = None # DataFrame. Columns {RA, DEC, Redshift, zbin}, rows = galaxies. bin_redshifs() adds column zbin

        # Mock box of random points
        self.Nrands = Nrands # number of random points to generate
        self.RandomsPS = None # DataFrame. Columns {x,y,z, vx,vy,vz, zbin}, rows = random points within a box defined by self.Lbox and self.zbox. zbin should be one of self.zbins
        self.RandomsRDZ = None # DataFrame. Columns {RA, DEC, Redshift, zbin}, rows = random points within a box defined by self.Lbox and self.zbox. zbin should be one of self.zbins

        # wtheta
        self.wtfout = wtfout # string. file name to save wtheta
        self.tbins = None # array. Theta [deg] at center of bin.
        self.tbin_edges = tbin_edges # array. Theta bin edges. if == None, set default below
        self.wtheta = None # DataFrame. Columns {theta bins}, rows wtheta, index zbin

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

    def getmock(self, Nstack=None, rtfout=None, zbin_width=None, Nrands=None, fow=None, nthreads=32, z4push=None, galplots=None):
        """
        Stacks Nstack^3 boxes together (around the origin) to create a bigger box.
            Nstack=0 => just moves origin to center of box in prep for push_box2catz.
        Pushes the box so the face is at comoving_distance(redshift = self.cat_zbox)
        Transforms to RA, DEC, Redshift and bins redshift using self.zbin_width.
        rtfout == string writes function runtimes to this file
               == None to skip timing fncs.
        fow == one of {None, 'wtheta', 'runtimes', 'all'}
                    => current file(s) will be renamed, new file(s) started

        Notes:

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
            self.report_times = OD([('mocknum', self.mocknum), ('nthreads',nthreads)]) # report_times dict for functions to report times
            self.report_times['fcol_width'] = 25 # set report_times file column width
            self.report_times['getmock'] = hf.time_code('start') #.TS. get code start time
            print('Function runtimes will be written to {}'.format(self.rtfout))
            print('getmock() started at {}'.format(datetime.datetime.now()))
        else:
            print('*** Function runtimes are NOT being tracked.\n\t Set MockBox.rtfout to track these.')

        #### The following are set on __init__, but can be changed here:
        if Nstack is not None:
            self.Nstack = Nstack

        if Nrands is not None:
            self.Nrands = Nrands

        if z4push is not None:
            self.z4push = z4push

        if zbin_width is not None:
            self.zbin_width = zbin_width

        # if statfout is not None:
        #     self.statfout = statfout
        #
        # if tbin_edges is not None:
        #     self.tbin_edges = tbin_edges
        # self.report_times['numtbins'] = len(self.tbin_edges)-1

        if galplots is not None:
            self.galplots = galplots
        ####

        if su.cosmo is None:
            su.load_cosmo() # loads default global cosmo object plus H0, Om0

        if fow is not None: # rename current files and start new ones
            self.ow_files(which=fow)
        ###


        # Get galaxy DF by populating DM mock using HODmodel
        self.cat_galtbl, self.cat_Lbox, self.cat_zbox = su.get_galtbl(getas='DF')
        if self.galplots: # plot a random subsample of galaxies
            mp.plot_galaxies(self.cat_galtbl, gal_frac=0.005, coords='xyz', title='Original Mock')

        # Stack boxes and push box face to appropriate redshift.
        self.transform_mock(box='PhaseSpace') # Sets self.PhaseSpace and self.RDZ

        # Get box of random points
        self.get_randoms()

        if self.rtfout is not None: # Currently this is not written to file
            self.report_times['getmock'] = hf.time_code(self.report_times['getmock'], unit='min') #.TE. replace start time with runtime in minutes
            print('\n\t{0}\ngetmock() ran for {1:.1f} minutes.\n'.format(datetime.datetime.now(), self.report_times['getmock']))

        return None




    def getmock_calcwtheta(self, tbin_edges=None, wtfout=None, Nstack=None, rtfout=None, zbin_width=None, Nrands=None, fow=None, nthreads=32, z4push=None, galplots=None):
        """
        Stacks Nstack^3 boxes together (around the origin) to create a bigger box.
            Nstack=0 => just moves origin to center of box in prep for push_box2catz.
        Pushes the box so the face is at comoving_distance(redshift = self.cat_zbox)
        Transforms to RA, DEC, Redshift and bins redshift using self.zbin_width.
        Calculates wtheta using tbins, nthreads and writes results to wtfout.
        Calculates runtime of each wtheta calculation and outputs info to rtfout.
        rtfout == string writes function runtimes to this file
               == None to skip timing fncs.
        fow == one of {None, 'wtheta', 'runtimes', 'all'}
                    => current file(s) will be renamed, new file(s) started

        Notes:

        Expect BAO at ~4.5 degrees for z~0.5 (see plots/zthetaBAO.png)
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
            self.report_times = OD([('mocknum', self.mocknum), ('nthreads',nthreads)]) # report_times dict for functions to report times
            self.report_times['fcol_width'] = 25 # set report_times file column width
            self.report_times['getmock_calcwtheta'] = hf.time_code('start') #.TS. get code start time
            print('Function runtimes will be written to {}'.format(self.rtfout))
            print('getmock_calcwtheta() started at {}'.format(datetime.datetime.now()))
        else:
            print('*** Function runtimes are NOT being tracked.\n\t Set MockBox.rtfout to track these.')

        #### The following are set on __init__, but can be changed here:
        if Nstack is not None:
            self.Nstack = Nstack

        if Nrands is not None:
            self.Nrands = Nrands

        if z4push is not None:
            self.z4push = z4push

        if zbin_width is not None:
            self.zbin_width = zbin_width

        if wtfout is not None:
            self.wtfout = wtfout

        if tbin_edges is not None:
            self.tbin_edges = tbin_edges
        self.report_times['numtbins'] = len(self.tbin_edges)-1

        if galplots is not None:
            self.galplots = galplots
        ####

        if su.cosmo is None:
            su.load_cosmo() # loads default global cosmo object plus H0, Om0

        if fow is not None: # rename current files and start new ones
            self.ow_files(which=fow)
        ###


        # Get galaxy DF by populating DM mock using HODmodel
        self.cat_galtbl, self.cat_Lbox, self.cat_zbox = su.get_galtbl(getas='DF')
        if self.galplots: # plot a random subsample of galaxies
            mp.plot_galaxies(self.cat_galtbl, gal_frac=0.005, coords='xyz', title='Original Mock')

        # Stack boxes and push box face to appropriate redshift.
        self.transform_mock(box='PhaseSpace') # Sets self.PhaseSpace and self.RDZ

        # Get box of random points and groupby redshift bins
        self.get_randoms()
        rgroups = self.RandomsRDZ.groupby('zbin', axis=0)

        # Calculate wtheta for each zbin and write to file
        zgroups = self.RDZ[['RA','DEC','zbin']].groupby('zbin', axis=0) # group by redshift bin, only need these columns
        for i, (zzz, rdz_z) in enumerate(zgroups):
            Randoms_z = rgroups.get_group(zzz) # get the randoms within this zbin

            # Several report_times entries are set in calc_write_wtheta and overwritten for each zzz.
            # Be sure to write report_times to file before the end of this for loop
            self.calc_write_wtheta(zzz, rdz_z, Randoms_z, nthreads=nthreads)

            if self.rtfout is not None: # Write report_times to file
                hf.write_report_times(self.report_times, self.rtfout)
                print('\tcalc_wtheta for zbin = {0} took {1:.1f} minutes'.format(zzz, self.report_times['calc_wtheta']))
                print('Results written to {0}. Calculation report_times written to {1}.\n'.format(self.wtfout, self.rtfout))
                # print('{0:15d} {1:15.1f} {2:15.1f} {3:15d} {4:15.4f} {5:15d}'.format(nthreads, zzz, ztime, len(rdz_z.index), dtm, len(tbins)), file=open(rtfout, 'a'))
                # print('\nwtheta calculation took {0:.1f} minutes with nthreads = {1}\n'.format(ztime, nthreads))

        if self.rtfout is not None: # Currently this is not written to file
            self.report_times['getmock_calcwtheta'] = hf.time_code(self.report_times['getmock_calcwtheta'], unit='min') #.TE. replace start time with runtime in minutes
            print('\n\t{0}\ndo_mock_wtheta.py ran for {1:.1f} minutes.\n'.format(datetime.datetime.now(), self.report_times['getmock_calcwtheta']))

        return None


    def ow_files(self, which='all'):
        # move current files so not overwritten

        if which == 'wtheta':
            mv_fout = hf.file_ow(self.wtfout)
        elif which == 'runtimes':
            mv_fout = hf.file_ow(self.rtfout)
        elif which == 'all':
            mv_fout = hf.file_ow(self.wtfout)
            mv_fout = hf.file_ow(self.rtfout)
        else:
            print('ow_files received invalid argument, which = {}'.format(which))

        return None


    def get_randoms(self, Nrands=None):
        """ Set self.RandomsPS with columns {x,y,z, vx,vy,vz, zbin}
            and self.RandomsRDZ with columns {RA, DEC, Redshift, zbin}, RA & DEC in degrees"""

        if Nrands is not None: # can be set on __init__ and/or changed here
            self.Nrands = Nrands
        if self.Nrands is None: # set a default
            self.Nrands = len(self.RDZ)
            print('MB.Nrands was not instantiated. Using Nrands = len(RDZ) = {}'.format(self.Nrands))

        if self.rtfout is not None:
            self.report_times['get_randoms'] = hf.time_code('start') #.TS. get code start time
            self.report_times['numrands'] = self.Nrands

        # create random points in box with side length self.Lbox, centered around origin
        ran_coords = np.random.random((self.Nrands,3))*self.Lbox - self.Lbox/2
        ran_vels = np.zeros((self.Nrands,3))
        # Set self.RandomsRDZ DF
        self.RandomsPS = pd.DataFrame(np.hstack([ran_coords,ran_vels]), columns=['x','y','z', 'vx','vy','vz'])
        self.push_box2catz(box='Randoms') # pushes the x-face to the catalog mock redshift

        # transform coords to RA, DEC, Redshift
        if self.rtfout is not None:
            self.report_times['get_ra_dec_z_Rands'] = hf.time_code('start')
        self.RandomsRDZ = hf.get_ra_dec_z(self.RandomsPS, usevel=False) # DF of {RA,DEC,Redshift}
        if self.rtfout is not None:
            self.report_times['get_ra_dec_z_Rands'] = hf.time_code(self.report_times['get_ra_dec_z_Rands'], unit='min') #.TE. replace start time with runtime in minutes

        # Find/set redshift bins
        self.bin_redshifs(box='Randoms', validate=False) # adds column 'zbin' to RandomsPS and RandomsRDZ
        if self.galplots:
            mp.plot_galaxies(self.RandomsPS, title="Galaxy Randoms at z_catalog")
            mp.plot_galaxies(self.RandomsPS, plotdim=2, title="2D: Galaxy Randoms at z_catalog")

        if self.rtfout is not None:
            self.report_times['get_randoms'] = hf.time_code(self.report_times['get_randoms']) # get function runtime [min]
            print('\tget_randoms() took {0:.1f} minutes'.format(self.report_times['get_randoms']))

        return None



    def calc_write_wtheta(self, zzz, rdz, randoms, wtfout=None, nthreads=32):
        """ zzz = Redshift at center of the mock or the mock slice
            rdz = DataFrame with min cols {RA, DEC} of galaxy positions
            randoms = DataFrame with min cols {RA, DEC} of random positions
                        in the same box or slice as rdz.
            wtfout = string, path of file to write wtheta to
        """

        if wtfout is not None: # this is set on __init__, but can be changed here
            self.wtfout = wtfout

        if self.rtfout is not None: # Save some info to report_times
            self.report_times['zbin'] = zzz
            self.report_times['numgals'] = len(rdz.index)
            print('\nCalculating wtheta for zbin = {0:1.2f}\n\t{1}\n'.format(zzz, datetime.datetime.now()))
            self.report_times['calc_wtheta'] = hf.time_code('start') #.TS. get code start time

        randoms_kwargs = { 'randbox':randoms }
        tbcens, wtheta, self.report_times = cw.calc_wtheta(rdz, MockBox=self, randoms_kwargs=randoms_kwargs, nthreads=nthreads)

        # Set self.tbins or check that it equals tbcens
        if self.tbins is None:
            self.tbins = tbcens
        else:
            errmsg = 'Theta bin centers from cw.calc_wtheta() do not match those previously set in MockBox.tbins.'
            np.testing.assert_allclose(tbcens, self.tbins, rtol=1e-5, err_msg=errmsg)
            # stops execution with error if theta bin centers to not match

        # Set or append to self.wtheta
        tmpdf = pd.DataFrame(data=wtheta, columns=[zzz], index=tbcens)
        tmpdf = tmpdf.T # so that theta bins are columns
        self.wtheta = tmpdf if self.wtheta is None else self.wtheta.append(tmpdf, ignore_index=False, sort=True)

        # Write current zbin wtheta info to file
        # Do this now so it is not lost if a future zbin calculation fails
        # cw.write_to_file(tbcens, wtheta, zzz, self.mocknum, self.wtfout)
        self.write_wtheta_to_file(tbcens, wtheta, zzz, len(rdz), len(randoms))

        if self.rtfout is not None:
            self.report_times['calc_wtheta'] = hf.time_code(self.report_times['calc_wtheta'], unit='min') #.TE. replace start time with runtime in minutes

        return None



    def write_wtheta_to_file(self, bcens, wtheta, zbin, Ngals, Nrands, fout=None):
        """ bcens = array of theta bin centers
            wtheta = array of wtheta values for each theta bin
            zbin = center of the redshift bin of galaxies used to calc wtheta
            fout = path (as string) of file to write or append to
                 = None uses self.wtfout

            Writes wtheta info to fout.
            Moves existing file if structure is not compatible.
        """
        # Setup
        if fout is not None: # this is set on __init__, but can be changed here
            self.wtfout = fout

        print('You should update this function (MockBox.write_wtheta_to_file) to print the proper PRECISION!')
        numbcens = len(bcens)
        # if you change this you must update several things in the rest of this function:
        extra_cols = np.array(['mock', 'zbin', 'Nstack', 'Ngals', 'Nrands'])

        # Check if file exists
        fpath = Path(self.wtfout)
        if fpath.is_file():
            rtol=1e-5
            # print('File {} exists. Checking compatibility...'.format(fout))
            # Check that bcens and extra_cols are the same as in current file.
            __, file_xcols, tbool = cw.get_tbins(self.wtfout, val_array=bcens, rtol=rtol)
            xbool = np.array_equal(file_xcols,extra_cols)
            if not (tbool and xbool): # columns don't match
                # Move the current file so we can start a new one.
                mv_fout = hf.file_ow(self.wtfout)
                print('*** thetabins' if not tbool else '***', 'extra_cols' if not xbool else '', 'not compatible with current file: {} ***'.format(self.wtfout))
                print('*** Moved existing file to {} so it is not overwritten. ***'.format(mv_fout))
            # else: # columns match
            #     print('Input data compatible (bcens rtol={0}) with existing file {1}.'.format(rtol, fout))

        # If fout has been moved (above) or never existed, create new file.
        if not fpath.is_file():
            # print('Writing new file {}'.format(fout))
            hdr = 'Columns labeled with floats contain bin centers (theta in degrees), others are extra info.\n'
            new_cols = np.stack([np.concatenate((extra_cols, bcens.astype(str)))])
            np.savetxt(self.wtfout, new_cols, fmt='%25.7s', header=hdr) # write header

        # Now fout exists, so append the data.
        print('Appending wtheta to {}'.format(self.wtfout))
        mlw = 50*(len(extra_cols) + len(bcens))
        dat_xcols = [self.mocknum, zbin, self.Nstack, Ngals, Nrands ]
        str_cols = np.array2string(np.append(dat_xcols, wtheta), \
                formatter={'float_kind':lambda x: "%25.15e" % x}, max_line_width=mlw)[1:-1]
        print(str_cols, file=open(self.wtfout, 'a'))

        return None





    # SDSS DR10 max photo zErr = 0.365106
    def bin_redshifs(self, box='RDZ', zbin_width=None, validate=False):
        """ Adds column 'zbin' to self.RDZ and self.PhaseSpace containing center of zbin the galaxy is in
                (zbin is one of self.zbins)
            Sets self.zbins and self.zbin_edges

            zbin_width = desired spacing between bins (result will not be exact)
        """

        # Get the right DataFrame
        if box=='RDZ':
            rtname = 'bin_redshifs'
            df = self.RDZ
            dfps = self.PhaseSpace
        elif box=='Randoms':
            rtname = 'bin_redshifs_Rands'
            df = self.RandomsRDZ
            dfps = self.RandomsPS
        else: # raise an error
            assert 0, 'bin_redshifs() received invaid argument: box = {}'.format(box)

        if self.rtfout is not None: # track function runtimes
                self.report_times[rtname] = hf.time_code('start') #.TS. get code start time

        if zbin_width is not None: # this is set on __init__, but can be changed here
            self.zbin_width = zbin_width

        print('\n*** You should fix redshift bins so you get consistent binning with different mocks. ***')

        try:
            self.zbin_edges
            self.zbins
            print('\nCheck this mock and make sure zBINS is a larger interval than zVALUES:')
            print('zBINS min, max = [{binmin},{binmax}]'.format(binmin=self.zbin_edges[0], binmax=self.zbin_edges[-1]))
            print('zVALUES min, max = [{zmin},{zmax}]'.format(zmin=self.RDZ.Redshift.min(), zmax=self.RDZ.Redshift.max()))
        except:
            tol = 0.03
            zmin = int(np.floor((self.RDZ.Redshift.min()-tol)*10))/10. # floor fnc (value-tol) in 1st decimal place
            zmin = max(0.0, zmin)
            zmax = int(np.ceil((self.RDZ.Redshift.max()+tol)*10))/10. # gives float to 1 decimal place
            num_binedges = int(np.ceil((zmax-zmin)/self.zbin_width))
            num_binedges = max(2, num_binedges)
            self.zbin_edges = np.linspace(zmin,zmax,num_binedges)
            self.zbins = np.round((self.zbin_edges[:-1]+self.zbin_edges[1:])/2, 2) # keep 2 decimal places
            print('\nMB.zbin_edges and MB.zbins have been set.')
            print('zBINS min, max = [{binmin},{binmax}]'.format(binmin=self.zbin_edges[0], binmax=self.zbin_edges[-1]))
            print('zVALUES min, max = [{zmin},{zmax}]'.format(zmin=self.RDZ.Redshift.min(), zmax=self.RDZ.Redshift.max()))

        # create zbin masks for rdz dataframe
            # add a column to rdz containing the zbin (self.zbins value) the galaxy resides in
        # given z, find which bin it's in and get the value of the bin center
        df['zbin'] = df['Redshift'].apply(hf.find_bin_center, **{"bin_edges":self.zbin_edges, "bin_centers":self.zbins})
        # set zbin in PhaseSpace df
        dfps['zbin'] = df['zbin']


        # make sure the operation worked as expected:
        if validate:
            for index, row in self.RDZ.iterrows():
                rdzbin = row.zbin
                truezbin = find_bin_center(row.Redshift, bin_edges=self.zbin_edges, bin_centers=self.zbins)
                assert(rdzbin==truezbin)


        if self.rtfout is not None: # track function runtimes
                self.report_times[rtname] = hf.time_code(self.report_times[rtname], unit='min') #.TE. replace start time with runtime in minutes

        return None





    def transform_mock(self, box='PhaseSpace'):
        """ Stack boxes around origin.
            Push face to z.
            Set self.PhaseSpace and self.RDZ

            All self.mock_* should be set prior to calling this.
                e.g. self.cat_galtbl, self.cat_Lbox, self.cat_zbox = su.get_galtbl(getas='DF')
         """

        # Stack boxes
        self.stack_boxes() # Sets self.PhaseSpace by stacking cat_galtbl, centered on origin.
        if self.galplots:
            mp.plot_galaxies(self.PhaseSpace, gal_frac=5e-4, coords='xyz', title='Boxes Stacked Around Origin')

        # Push box face to redshift = self.cat_zbox
        self.push_box2catz(box=box)
        if self.galplots:
            mp.plot_galaxies(self.PhaseSpace, gal_frac=5e-4, coords='xyz', title='Boxes Stacked and Pushed to Redshift = {}'.format(self.zbox))

        # Transform coordinates to RA, DEC, Redshift
        print('\nConverting to RA, DEC, Redshift. ...')
        if self.rtfout is not None:
            self.report_times['get_ra_dec_z'] = hf.time_code('start')
        self.RDZ = hf.get_ra_dec_z(self.PhaseSpace, usevel=True)
        if self.rtfout is not None:
            self.report_times['get_ra_dec_z'] = hf.time_code(self.report_times['get_ra_dec_z'], unit='min') #.TE. replace start time with runtime in minutes

        # Bin redshifts, add column to self.RDZ and self.PhaseSpace
        self.bin_redshifs(box='RDZ', validate=False)

        # Plot galaxy distribution
        if self.galplots:
            # plot colored by zbins
            mp.plot_galaxies(pd.concat([self.PhaseSpace,self.RDZ['zbin']],axis=1), title="Final Galaxies box colored by zbin")
            mp.plot_galaxies(pd.concat([self.PhaseSpace,self.RDZ['zbin']],axis=1), plotdim=2, title="2D: Final Galaxies box colored by zbin")

        return None



    def stack_boxes(self):
        """
        uses galdf = self.cat_galtbl: DataFrame representing a mock box,
                with minimum columns {'x','y','z', 'vx','vy','vz'},
                (assumes coordinates {x,y,z} in [0,Lbox]).
        Stacks self.Nstack^3 (must be an EVEN integer) galdf mock boxes together around the origin.
        Sets self.Lbox.
        Sets self.PhaseSpace: DataFrame of galaxies with the origin at the center of the box and boxsize=self.Lbox.
                columns {x,y,z} generated periodically, {vx,vy,vz, all others} same as original galaxies.
        """
        print('\nStacking {}^3 boxes. ...'.format(self.Nstack))

        # Setup
        if self.rtfout is not None: # track function runtimes
            self.report_times['stack_boxes'] = hf.time_code('start') #.TS. get code start time

        self.Lbox = self.cat_Lbox*self.Nstack

        # check quadrant
        galdf = self.cat_galtbl
        print('*** Warning: MockBox.stack_boxes assumes the original box is strictly in the 1st quadrant with the origin at the corner. ***')
        print('\tVerifying that original coordinates are non-negative (this fnc does not check for a positive offset from 0).')
        assert sum(pd.concat([galdf.x<0, galdf.y<0, galdf.z<0])) == 0, \
                'Quadrant test failed: Some MockBox.cat_galtbl coordinates are not in the 1st quadrant as required.'

        # Stack boxes
        if self.Nstack == 0: # Just move the origin to the center of the box.
            print('Moving origin to box center...')
            self.Lbox = self.cat_Lbox
            L2 = self.Lbox/2.
            newgals = galdf.copy(deep=True)
            newgals.x = newgals.x-L2; newgals.y = newgals.y-L2; newgals.z = newgals.z-L2

        else:
            newgals = self.stack_boxes_doit()
            # index of original galaxy is retained => each index is repeated Nstack^3 times.

        self.PhaseSpace = newgals

        if self.rtfout is not None: # track function runtimes
            self.report_times['stack_boxes'] = hf.time_code(self.report_times['stack_boxes'], unit='min') #.TE. replace start time with runtime in minutes
            print('\n\t{0}\nstack_boxes() ran for {1:.1f} minutes.\n'.format(datetime.datetime.now(), self.report_times['stack_boxes']))

        return None


    def stack_boxes_doit(self):
        """ Uses df = self.cat_galtbl
            Generates new 'xyz' coordinates by stacking Nstack**3 boxes.
            Copy all other columns.
            Returns a DataFrame with Nstack^3 galaxies (rows)
        """
        df = self.cat_galtbl
        Nstack = self.Nstack

        N2 = int(Nstack/2)
        assert Nstack%2 == 0, 'Nstack should be an even integer.'

        dfcols = list(df.columns.values)
        extra_vals = df[[c for c in dfcols if c not in ['x','y','z']]] # get df of just extra values

        idx_omag = np.floor(np.log10(np.max(df.index.values))) # order of magnitude of max df.index.value
        nblist = [] # list to hold DataFrames for each new box
        nLin = np.linspace(-N2,N2-1,Nstack)*self.cat_Lbox  # array of ints -N2 to N2-1, defines deltas in each direction
        deltax,deltay,deltaz = np.meshgrid(nLin,nLin,nLin) # holds deltas for each new box
        for b in range(Nstack**3): # for each new box
            boxdf = extra_vals.copy(deep=True)
            dx,dy,dz = deltax.flat[b], deltay.flat[b], deltaz.flat[b]
            boxdf['x'] = df.x + dx
            boxdf['y'] = df.y + dy
            boxdf['z'] = df.z + dz

            idx_offset = np.int32(b*10**(idx_omag+1)) # add this to the indices to ensure unique indices for each galaxy
            boxdf.index = boxdf.index + idx_offset

            nblist.append(boxdf)

        newdf = pd.concat(nblist, ignore_index=False)
        # index of original galaxy is retained => each index is repeated Nstack^3 times.

        return newdf


    def stack_boxes_gen_new_rows_old(self, gal):#, Nstack=None, ogLbox=None):
        """ This gets applied to each row of the original dataframe.
            gal = single galaxy data as series (single row from original df)
            Generates new 'xyz' coordinates by stacking self.Nstack**3 boxes.
            Copy all other columns.
            Returns a DataFrame with self.Nstack**3 galaxies (rows)
        """
        N2 = int(self.Nstack/2)
        assert self.Nstack%2 == 0, 'Nstack should be an even integer.'

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



    def stack_boxes_old(self):
        """
        uses galdf = self.cat_galtbl: DataFrame representing a mock box,
                with minimum columns {'x','y','z', 'vx','vy','vz'},
                (assumes coordinates {x,y,z} in [0,Lbox]).
        Stacks self.Nstack^3 (must be an EVEN integer) galdf mock boxes together around the origin.
        Sets self.Lbox.
        Sets self.PhaseSpace: DataFrame of galaxies with the origin at the center of the box and boxsize=self.Lbox.
                columns {x,y,z} generated periodically, {vx,vy,vz, all others} same as original galaxies.
        """
        # Setup
        if self.rtfout is not None: # track function runtimes
            self.report_times['stack_boxes'] = hf.time_code('start') #.TS. get code start time

        self.Lbox = self.cat_Lbox*self.Nstack
        print('\nStacking {}^3 boxes. ...'.format(self.Nstack))
        print('*** Warning: MockBox.stack_boxes assumes the original box is strictly in the 1st quadrant with the origin at the corner. ***')

        if self.Nstack == 0: # Just move the origin to the center of the box.
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
                nglist.append(self.stack_boxes_gen_new_rows(gal))
            # Create new df from dfs in nglist
            newgals = pd.concat(nglist, ignore_index=True)

        self.PhaseSpace = newgals

        if self.rtfout is not None: # track function runtimes
            self.report_times['stack_boxes'] = hf.time_code(self.report_times['stack_boxes'], unit='min') #.TE. replace start time with runtime in minutes
            print('\n\t{0}\nstack_boxes() ran for {1:.1f} minutes.\n'.format(datetime.datetime.now(), self.report_times['stack_boxes']))

        return None



    # cosmo.comoving_distance(zred) returns "Comoving distance in Mpc to each input redshift.""
    # from help(cosmo): Dimensionless Hubble constant: h = H_0 / 100 [km/sec/Mpc]
    def push_box2catz(self, box='PhaseSpace'):
        """ self.PhaseSpace (or self.RandomsRDZ depending on box=) must be set before calling this function.
                ASSUMES the box is currently centered at the origin!
            Moves the box coordinates so the x-FACE is at comoving_distance(self.cat_zbox).
                (Only 'x' column is changed.)
            Sets self.zbox = self.cat_zbox if box=='PhaseSpace'
        """
        # global cosmo
        if su.cosmo is None:
            su.load_cosmo()

        # Get the redshift for the push
        try:
            z4push = self.cat_zbox if self.z4push=='cat' else np.float32(self.z4push)
        except ValueError as ve:
            print("MockBox.z4push should be either 'cat' or a float.\nReceived unrecognized option z4push = {}".format(self.z4push))
            raise ve

        # Get the right DataFrame
        if box=='PhaseSpace':
            rtname = 'push_box2catz'
            self.zbox = z4push
            df = self.PhaseSpace
        elif box=='Randoms':
            rtname = 'push_box2catz_Rands'
            df = self.RandomsPS
        else: # raise an error
            assert 0, 'push_box2catz() received invaid argument: box = {}'.format(box)
        assert df is not None, 'You must instantiate MB.{box} before calling push_box2catz().'.format(box=box)

        print('\nPushing the {b} box out to box x-face redshift = {z:1.2f} ...'.format(b=box, z=z4push))
        if self.rtfout is not None: # track function runtimes
            self.report_times[rtname] = hf.time_code('start') #.TS. get code start time

        # Shift x so coordinates are strictly positive (i.e. move face to x=0)
        # Then push the face to comoving_distance(self.cat_zbox)
        deltax = self.Lbox/2. + (su.cosmo.comoving_distance(z4push).value)*su.cosmo.h # Mpc/h
        df.x = df.x + deltax

        if self.rtfout is not None: # track function runtimes
            self.report_times[rtname] = hf.time_code(self.report_times[rtname], unit='min') #.TE. replace start time with runtime in minute

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
