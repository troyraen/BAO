import numpy as np
from pathlib import Path

from halotools.mock_observables import mock_survey
# import Corrfunc
from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks
from Corrfunc.utils import convert_3d_counts_to_cf

import setup_mock as sm


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
        fout = path (as string) of file to write or append to

        Writes file fout:
            First line is bcens values, then 'zbin' and 'mock' number
            All following lines are wtheta values for corresponding bcen, then 'zbin' and 'mock'
    """
    print('\nYou should update this function (calc_wtheta.write_to_file) to print the proper PRECISION!\n')
    extra_cols = np.array(['zbin', 'mock'])

    # check if file exists
    fpath = Path(fout)
    if fpath.is_file(): # if it does, check column compatibility
        print('File {} exists. Checking compatibility...'.format(fout))
        # check that bcens are the same
        file_cols = np.genfromtxt(fout, max_rows=1, dtype=str)
        numbcens = len(bcens)
        file_bcens = file_cols[:numbcens].astype(np.double)
        rtol=1e-5
        np.testing.assert_allclose(file_bcens, bcens, rtol=rtol, \
            err_msg='\nbcens != bcens (first line) in {}.\nwtheta not written to file.\n'.format(fout))
        # check that extra_cols are the same
        file_xcols = file_cols[numbcens:]
        assert np.array_equal(file_xcols,extra_cols), "\nextra columns don't match (first line) in {}.\nwtheta not written to file.\n".format(fout)
        # if they are, then append to file (below)
        print('Input data compatible with existing file (bcens rtol={}). Appending wtheta.\n'.format(rtol))

    else: # else create new file and write header
        print('Writing new file {}'.format(fout))
        hdr = 'First row contains bin centers (theta in degrees), then extra info. All other rows contain wtheta for that bin, then extra info.\n'
        new_cols = np.stack([np.concatenate((bcens.astype(str), extra_cols))])
        np.savetxt(fout, new_cols, fmt='%25.7s', header=hdr) # write header

    # append data to the file
    srtln = 50*(len(wtheta)+len(extra_cols))
    dat_cols = np.append(wtheta, [zbin,mocknum]) # add extra column data to wtheta
    str_cols = np.array2string(dat_cols, formatter={'float_kind':lambda x: "%25.15e" % x}, max_line_width=srtln)[1:-1]
    print(str_cols, file=open(fout, 'a'))




def load_from_file(fin):
    # returns pandas dataframe with bin centers as column headers
    # wtheta in rows (1 mock per row)
    return pd.read_csv(fin, delim_whitespace=True, comment='#')






def calc_wtheta(galaxy_df, bins, randoms_kwargs, nthreads=48):
    """galaxy_df = DataFrame including (at least) columns 'RA' and 'DEC'
        bins = array of theta bin edges in degrees
        randoms_kwargs (for get_randoms) can include: {Nran=, boxsize=, push_to_z=, cosmo=}
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
    RA, DEC = galaxy_df.RA, galaxy_df.DEC
    # CHECK PLOT TO MAKE SURE CORRD TRANSFORM IS AS EXPECTED
    DD_counts = DDtheta_mocks(autocorr, nthreads, bins, RA, DEC)

    # random random
    autocorr=1
    if 'Nran' not in randoms_kwargs:
        randoms_kwargs['Nran'] = len(RA)*10
    print('\nGetting randoms with')
    print(randoms_kwargs)
    rand_RA, rand_DEC = get_randoms(**randoms_kwargs)
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




def get_randoms(Nran=10**5, boxsize=1000, push_to_z=None, cosmo=None):
    """Returns random [RA, DEC] in degrees"""
    ran_coords = np.random.random((Nran,3))*boxsize
    ran_vels = np.zeros((Nran,3))
    ps_coords = np.hstack([ran_coords,ran_vels])
    if push_to_z is not None:
        ps_coords = sm.push_box2z(ps_coords, push_to_z, boxsize, cosmo=cosmo) # returns original ndarray with 1st column shifted
    ran_ra, ran_dec, ran_z = get_ra_dec_z(ps_coords, cosmo=cosmo, usevel=True)

# using Duncan's function:
# https://halotools.readthedocs.io/en/latest/_modules/halotools/mock_observables/mock_survey.html
    # ran_ra, ran_dec, ran_z = mock_survey.ra_dec_z(ran_coords, ran_vels)
    # ran_ra = np.degrees(ran_ra) # [0, 90] degrees
    # ran_dec = np.degrees(ran_dec) # [-90, 0] degrees

    return [ran_ra, ran_dec]
