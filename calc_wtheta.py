import numpy as np
from pathlib import Path

from halotools.mock_observables import mock_survey
# import Corrfunc
from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks
from Corrfunc.utils import convert_3d_counts_to_cf



def get_wtheta(halocat, HODmodel, bins, repop=True, fout=None):
    """Takes a halo catalog and HOD model (repopulates the mock if repop==True).
        bins = array of bin edges in degrees
        Pass fout = 'file_path' to write wtheta to a file (will append if file exists)
    Returns statistics: wtheta in given bins
    """

    if repop: HODmodel.mock.populate()
    galtbl = HODmodel.mock.galaxy_table
    boxsize = halocat.Lbox[0]
    # redshift = halocat.redshift
    bcens, wtheta = calc_wtheta(galtbl, bins, boxsize=boxsize)

    if type(fout) == str:
        write_to_file(bcens, wtheta, fout)

    return [bcens, wtheta]





#### HELPER FUNCTIONS ####

def write_to_file(bcens, wtheta, fout):
    # check if file exists
    fpath = Path(fout)
    if fpath.is_file():
        print('File {} exists. Checking compatibility...'.format(fout))
    # if it does, check that bcens are the same
        bcens0 = np.genfromtxt(fout, max_rows=1)
        np.testing.assert_allclose(bcens0, bcens, \
            err_msg='\nbcens != bcens (first line) in {}.\nwtheta not written to file.\n'.format(fout))
        print('bcens compatible with existing file. Appending wtheta...')

        srtln = 50*len(wtheta)
        wtstr = np.array2string(wtheta, formatter={'float_kind':lambda x: "%25.15e" % x}, max_line_width=srtln)[1:-1]
        print(wtstr, file=open(fout, 'a'))
        #
        # f = open(fout, 'ab')
        # # data = np.stack([wtheta])
        # np.savetxt(f, wtheta, fmt='%25.15e')
    # , then append to file
    # else:
    hdr = 'First row contains bin centers. All other rows contain wtheta for that bin.\n'
    data = np.stack([bcens, wtheta])
    np.savetxt(fout, data, fmt='%25.15e', header=hdr)
    # np.savetxt(fout, data, fmt=['%20.9f', '%20.9e'], header=hdr)
    print('\nYou should update this function (calc_wtheta.write_to_file) to print the proper PRECISION!\n')


# def load_from_file(fin):


def get_ra_dec_z(galaxy_table):
    coords = np.vstack([galaxy_table['x'], galaxy_table['y'], galaxy_table['z']]).T # check these units, mock_survey.ra_dec_z expects Mpc/h
    vels = np.vstack([galaxy_table['vx'], galaxy_table['vy'], galaxy_table['vz']]).T # mock_survey.ra_dec_z expects km/s

    ra, dec, z = mock_survey.ra_dec_z(coords, vels) # returns ra, dec in radians
    ra = np.degrees(ra) # [0, 90] degrees
    dec = np.degrees(dec) # [-90, 0] degrees

    return [ra, dec, z]



def calc_wtheta(galaxy_table, bins, nthreads=48, boxsize=1000):
    """bins = array of bin edges in degrees
    Returns [bcens, wtheta]
    """

    # Calc pair counts---
    # bins = np.linspace(binrange[0], binrange[1], nbins + 1) # theta values [degrees] for bin-edges

    # DDtheta_mocks expects ra in [0.0, 360.0] and dec in [-90.0, 90.0] degrees
    # DDtheta_mocks returns numpy structured array containing
    # [thetamin, thetamax, thetaavg, npairs, weightavg] for each angular bin

    # gal gal
    autocorr=1
    RA, DEC, __ = get_ra_dec_z(galaxy_table)
    # CHECK PLOT TO MAKE SURE CORRD TRANSFORM IS AS EXPECTED
    DD_counts = DDtheta_mocks(autocorr, nthreads, bins, RA, DEC)

    # random random
    autocorr=1
    Nran = len(RA)*10
    rand_RA, rand_DEC = get_randoms(Nran=Nran, boxsize=boxsize)
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




def get_randoms(Nran=10**5, boxsize=1000):
    """Returns random [RA, DEC] in degrees, as float32 to match output from get_ra_dec_z()"""
    ran_coords = np.random.random((Nran,3))*boxsize
    ran_vels = np.zeros((Nran,3))

# https://halotools.readthedocs.io/en/latest/_modules/halotools/mock_observables/mock_survey.html
    ran_ra, ran_dec, ran_z = mock_survey.ra_dec_z(ran_coords, ran_vels)
    ran_ra = (np.degrees(ran_ra)).astype('float32') # [0, 90] degrees
    ran_dec = (np.degrees(ran_dec)).astype('float32') # [-90, 0] degrees

    return [ran_ra, ran_dec]
