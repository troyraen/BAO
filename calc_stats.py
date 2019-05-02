import numpy as np
import pandas as pd
from pathlib import Path
from astropy.table import Table
# import os
import datetime

from halotools.mock_observables import mock_survey
# import Corrfunc
from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks
from Corrfunc.utils import convert_3d_counts_to_cf
from Corrfunc.theory.wp import wp
from Corrfunc.theory.xi import xi

import setup as su
import helper_fncs as hf
import myplots as mp



def calc_wp(MockBox, nthreads=24):

    boxsize = MockBox.Lbox
    bins = MockBox.rbin_edges
    pimax = MockBox.pimax
    rt = MockBox.report_times
    print('\nCalculating wp(rp)\n\t{1}\n'.format(datetime.datetime.now()))

    # MockBox has x-face pushed to catalog redshift.
    # Must transform this to z direction.
    X, Y, Z = MockBox.PhaseSpace.y, MockBox.PhaseSpace.z, MockBox.PhaseSpace.x

    rt['calc_wp'] = hf.time_code('start')
    results = wp(boxsize, pimax, nthreads, bins, X, Y, Z)
    rt['calc_wp'] = hf.time_code(rt['calc_wp'], unit='min')

    bcens = np.around((bins[:-1]+bins[1:])/2, decimals=5)
    return [bcens, results['wp']]


def calc_xi(MockBox, nthreads=24):

    boxsize = MockBox.Lbox
    bins = MockBox.rbin_edges
    rt = MockBox.report_times
    print('\nCalculating xi\n\t{1}\n'.format(datetime.datetime.now()))

    # MockBox has x-face pushed to catalog redshift.
    # Transform this to z direction for consistency with calc_wp().
    X, Y, Z = MockBox.PhaseSpace.y, MockBox.PhaseSpace.z, MockBox.PhaseSpace.x

    rt['calc_xi'] = hf.time_code('start')
    results = xi(boxsize, nthreads, bins, X, Y, Z)
    rt['calc_xi'] = hf.time_code(rt['calc_xi'], unit='min')

    bcens = np.around((bins[:-1]+bins[1:])/2, decimals=5)
    return [bcens, results['xi']]



def calc_wtheta(galaxy_df, randoms_df, MockBox=None, nthreads=24):
    """ Pass a MockBox instance to use the following variables:
            tbins = MockBox.tbin_edges
            report_times = MockBox.report_times
        galaxy_df = galaxies DataFrame including (at least) columns 'RA' and 'DEC'
        randoms_df = randoms DataFrame including (at least) columns 'RA' and 'DEC'

    Returns [theta bin centers, wtheta in each bin, report_times]
    """

    rt = MockBox.report_times
    print('\nCalculating wtheta for zbin = {0:1.2f}\n\t{1}\n'.format(zzz, datetime.datetime.now()))
    rt['calc_wtheta'] = hf.time_code('start')

    RA, DEC = np.asarray(galaxy_df.RA), np.asarray(galaxy_df.DEC)
    N = len(RA)
    rand_RA, rand_DEC = np.asarray(randoms_df.RA), np.asarray(randoms_df.DEC)
    rand_N = len(rand_RA)

    tbins = MockBox.tbin_edges
    if tbins is None:
        tbins = np.logspace(np.log10(1.0), np.log10(10.0), 25)
        print('*** No theta bin edges provided to calc_wtheta.\n\t Using default tbins = {}'.format(tbins))

    try: # make sure nthreads was passed through correctly
        assert nthreads == rt['nthreads']
    except:
        print('\n*** WARNING: calc_wtheta() using nthreads={cwn}\n\tThis is different than set in MB.report_times.nthreads={rtn}\n'.format(cwn=nthreads, rtn=rt['nthreads']))


    #--- Calc pair counts---

    # DDtheta_mocks expects ra in [0.0, 360.0] and dec in [-90.0, 90.0] degrees
    # DDtheta_mocks returns numpy structured array containing
    #   [thetamin, thetamax, thetaavg, npairs, weightavg] for each angular bin

    # gal gal
    autocorr=1
    print('Calculating DD_counts...')
    rt['wt_galgal_counts'] = hf.time_code('start')
    DD_counts = DDtheta_mocks(autocorr, nthreads, tbins, RA, DEC)
    rt['wt_galgal_counts'] = hf.time_code(rt['wt_galgal_counts'], unit='min')


    # random random
    autocorr=1
    print('Calculating RR_counts...')
    rt['wt_randrand_counts'] = hf.time_code('start')
    RR_counts = DDtheta_mocks(autocorr, nthreads, tbins, rand_RA, rand_DEC)
    rt['wt_randrand_counts'] = hf.time_code(rt['wt_randrand_counts'], unit='min')


    # gal random
    autocorr=0
    print('Calculating DR_counts...')
    rt['wt_galrand_counts'] = hf.time_code('start')
    DR_counts = DDtheta_mocks(autocorr, nthreads, tbins, RA, DEC, RA2=rand_RA, DEC2=rand_DEC)
    rt['wt_galrand_counts'] = hf.time_code(rt['wt_galrand_counts'], unit='min')
    #---

    # calc w(theta)
    print('Converting counts to correlation...')
    rt['wt_counts_to_cf'] = hf.time_code('start')
    wtheta = convert_3d_counts_to_cf(N, N, rand_N, rand_N, DD_counts, DR_counts, DR_counts, RR_counts)
    rt['wt_counts_to_cf'] = hf.time_code(rt['wt_counts_to_cf'], unit='min')

    rt['calc_wtheta'] = hf.time_code(rt['calc_wtheta'], unit='min')

    bcens = np.around((tbins[:-1]+tbins[1:])/2, decimals=5)
    return [bcens, wtheta]
