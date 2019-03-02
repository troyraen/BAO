import numpy as np
import pandas as pd
from pathlib import Path
from astropy.table import Table
import os
import datetime

from halotools.mock_observables import mock_survey
# import Corrfunc
from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks
from Corrfunc.utils import convert_3d_counts_to_cf

import setup as su
import helper_fncs as hf
import myplots as mp


#### MAIN FUNCTION ####

def get_wtheta(halocat, HODmodel, bins, repop=True, fout=None):
    """Takes a halo catalog and HOD model (repopulates the mock if repop==True).
        bins = array of bin edges in degrees
        Pass fout = 'file_path' to write wtheta to a file (will append if file exists)
    Returns statistics: wtheta in given bins
    """

    print("\nDO NOT USE THIS FUNCTION (calc_wtheta.get_wtheta) BEFORE FIXING IT!")

    if repop: HODmodel.mock.populate()
    galtbl = HODmodel.mock.galaxy_table
    boxsize = halocat.Lbox[0]
    # redshift = halocat.redshift
    bcens, wtheta = calc_wtheta(galtbl, bins, boxsize=boxsize)

    try:
        write_to_file(bcens, wtheta, fout)
    except TypeError as e:
        print('\nNo file path given.\nbcens, wtheta not written to file.')
        print('bcens = {0}\nwtheta = {1}\n'.format(bcens, wtheta))
    except AssertionError as e:
        print('bcens = {0}\nwtheta = {1}\n'.format(bcens, wtheta))
        raise e

    finally:
        return [bcens, wtheta]





#### HELPER FUNCTIONS ####

def write_to_file(bcens, wtheta, zbin, mocknum, fout):
    """ bcens = array of theta bin centers
        wtheta = array of wtheta values for each theta bin
        zbin = center of the redshift bin of galaxies used for wtheta
        mocknum = float to use as mock ID number
        fout = path (as string) of file to write or append to

        Writes wtheta info to fout.
        Moves existing file if structure is not compatible.
    """
    # Setup
    print('\nYou should update this function (calc_wtheta.write_to_file) to print the proper PRECISION!\n')
    extra_cols = np.array(['zbin', 'mock'])
    numbcens = len(bcens)
    fpath = Path(fout)

    # Check if file exists
    if fpath.is_file():
        rtol=1e-5
        # print('File {} exists. Checking compatibility...'.format(fout))
        # Check that bcens and extra_cols are the same as in current file.
        __, file_xcols, tbool = get_tbins(fout, val_array=bcens, rtol=rtol)
        xbool = np.array_equal(file_xcols,extra_cols)
        if not (tbool and xbool): # columns don't match
            # Move the current file so we can start a new one.
            dtm = datetime.datetime.now() # get date and time to use in file name
            mv_fout = fout[:-4] + dtm.strftime("_ow_%m%d%y_%H%M") + fout[-4:] # assumes fout file extension is 3 letters
            os.rename(fout, mv_fout)
            print('*** thetabins' if not tbool else '***', 'extra_cols' if not xbool else '', 'not compatible with current file: {} ***'.format(fout))
            print('*** Moved existing file to {} so it is not overwritten. ***'.format(mv_fout))
        else: # columns match
            print('Input data compatible (bcens rtol={0}) with existing file {1}.'.format(rtol, fout))

    # If fout has been moved (above) or never existed, create new file.
    if not fpath.is_file():
        print('Writing new file {}'.format(fout))
        hdr = 'Columns labeled with floats contain bin centers (theta in degrees), others are extra info.\n'
        new_cols = np.stack([np.concatenate((bcens.astype(str), extra_cols))])
        np.savetxt(fout, new_cols, fmt='%25.7s', header=hdr) # write header

    # Now fout exists, so append the data.
    print('Appending wtheta to {}'.format(fout))
    mlw = 50*(len(bcens)+len(extra_cols))
    dat_cols = np.append(wtheta, [zbin,mocknum]) # add extra column data to wtheta
    str_cols = np.array2string(dat_cols, formatter={'float_kind':lambda x: "%25.15e" % x}, max_line_width=mlw)[1:-1]
    print(str_cols, file=open(fout, 'a'))




def load_from_file(fin):
    # returns pandas dataframe with bin centers as column headers
    # wtheta in rows (1 mock per row)
    return pd.read_csv(fin, delim_whitespace=True, comment='#')


def get_tbins(wdat, val_array=None, rtol=1e-5):
    """
    wdat = DataFrame or a file path as a string
    Assumes all column names that can be converted to floats are theta bins.
    bincols = list of theta bin column names as strings,
    ocols = list of other column names as strings
    If val_array is None
        Returns [ bincols, ocols ]
    Else
        Returns [bincols, ocols, boolean = (wdat tbins == val_array)]
    """

    wdf = cw.load_from_file(wdat) if (type(wdat) == str) else wdat
    allcols = list(wdf.columns.values)
    bincols = []
    ocols = []
    for c in allcols:
        try:
            float(c)
            bincols.append(c)
        except:
            ocols.append(c)

    if val_array is None:
        return [ bincols, ocols ]

    else:
        barr = np.asarray(bincols, dtype=np.double) # convert to array
        try:
            val_array = np.asarray(val_array, dtype=np.double) # try to convert
            np.testing.assert_allclose(barr, val_array, rtol=rtol )
        except TypeError:
            print('\nTypeError: val_array must be of type that can be converted to np.asarray(type=np.double)')
            return [bincols, ocols, False]
        except: # the bincols arrays are not equal
            print('\nbincols from input != val_array.')
            return [bincols, ocols, False]
        else:
            return [bincols, ocols, True]




def calc_wtheta(galaxy_df, bins, randoms_kwargs, nthreads=48):
    """galaxy_df = DataFrame including (at least) columns 'RA' and 'DEC'
        bins = array of theta bin edges in degrees
        randoms_kwargs (for get_randoms) can include: {Nran=, boxsize=, push_to_z=}
    Returns [theta bin centers, wtheta]
    """

    # Calc pair counts---
    # bins = np.linspace(binrange[0], binrange[1], nbins + 1) # theta values [degrees] for bin-edges

    # DDtheta_mocks expects ra in [0.0, 360.0] and dec in [-90.0, 90.0] degrees
    # DDtheta_mocks returns numpy structured array containing
    # [thetamin, thetamax, thetaavg, npairs, weightavg] for each angular bin

    # gal gal
    autocorr=1
    # RA, DEC, __ = get_ra_dec_z(galaxy_table)
    RA, DEC = np.asarray(galaxy_df.RA), np.asarray(galaxy_df.DEC)
    # CHECK PLOT TO MAKE SURE CORRD TRANSFORM IS AS EXPECTED
    DD_counts = DDtheta_mocks(autocorr, nthreads, bins, RA, DEC)

    # random random
    autocorr=1
    if 'Nran' not in randoms_kwargs:
        randoms_kwargs['Nran'] = len(RA)*10
    print('Getting randoms with {}'.format(randoms_kwargs))
    rand_RA, rand_DEC = get_randoms(**randoms_kwargs)
    print('Calculating...')
    RR_counts = DDtheta_mocks(autocorr, nthreads, bins, rand_RA, rand_DEC)

    # gal random
    autocorr=0
    DR_counts = DDtheta_mocks(autocorr, nthreads, bins, RA, DEC, RA2=rand_RA, DEC2=rand_DEC)
    #---

    # calc w(theta)
    N = len(RA)
    rand_N = len(rand_RA)
    wtheta = convert_3d_counts_to_cf(N, N, rand_N, rand_N, DD_counts, DR_counts, DR_counts, RR_counts)

    bcens = (bins[:-1]+bins[1:])/2
    return [bcens, wtheta]




def get_randoms(Nran=10**5, boxsize=1000, push_to_z=None, viewgals=False):
    """Returns random [RA, DEC] in degrees"""
    # create random points in box with side length boxsize, centered around origin
    ran_coords = np.random.random((Nran,3))*boxsize - boxsize/2

    ran_vels = np.zeros((Nran,3))
    ps_coords = np.hstack([ran_coords,ran_vels])
    # convert to DF to work with my functions
    pscdf = pd.DataFrame(ps_coords, columns=['x','y','z', 'vx','vy','vz'])
    if push_to_z is not None:
        pscdf = su.push_box2z(pscdf, push_to_z, boxsize) # returns original ndarray with 1st column shifted
        if viewgals:
            # plot to check coords
            # ngtbl = Table(ps_coords, names=['x','y','z', 'vx','vy','vz'])
            mp.plot_galaxies(pscdf, gal_frac=5e-5, coords='xyz', title="Galaxy Randoms")

    pscdf = hf.get_ra_dec_z(pscdf, usevel=True) # returns a DataFrame
    ran_ra, ran_dec = np.asarray(pscdf.RA), np.asarray(pscdf.DEC)
    # print('numgals = {0}, len(ran_ra) = {1}, len(ran_dec) = {2}'.format(Nran, len(ran_ra), len(ran_dec)))
    # print('Sampling ran_ra, ran_dec: {}'.format(ran_ra[:10], ran_dec[:10]))
# using Duncan's function:
# https://halotools.readthedocs.io/en/latest/_modules/halotools/mock_observables/mock_survey.html
    # ran_ra, ran_dec, ran_z = mock_survey.ra_dec_z(ran_coords, ran_vels)
    # ran_ra = np.degrees(ran_ra) # [0, 90] degrees
    # ran_dec = np.degrees(ran_dec) # [-90, 0] degrees

    return [ran_ra, ran_dec]
