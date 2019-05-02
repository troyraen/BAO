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

import setup as su
import helper_fncs as hf
import myplots as mp






def calc_wtheta(galaxy_df, randoms_df, MockBox=None, nthreads=32):
    """ Pass a MockBox instance to use the following variables:
            tbins = MockBox.tbin_edges
            report_times = MockBox.report_times
        galaxy_df = galaxies DataFrame including (at least) columns 'RA' and 'DEC'
        randoms_df = randoms DataFrame including (at least) columns 'RA' and 'DEC'

    Returns [theta bin centers, wtheta in each bin, report_times]
    """

    RA, DEC = np.asarray(galaxy_df.RA), np.asarray(galaxy_df.DEC)
    N = len(RA)
    rand_RA, rand_DEC = np.asarray(randoms_df.RA), np.asarray(randoms_df.DEC)
    rand_N = len(rand_RA)

    tbins = MockBox.tbin_edges
    if tbins is None:
        tbins = np.logspace(np.log10(1.0), np.log10(10.0), 25)
        print('*** No theta bin edges provided to calc_wtheta.\n\t Using default tbins = {}'.format(tbins))

    rt = MockBox.report_times
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
    rt['wt_galgal_counts'] = hf.time_code('start') #.TS. get code start time
    DD_counts = DDtheta_mocks(autocorr, nthreads, tbins, RA, DEC) #.TC.
    rt['wt_galgal_counts'] = hf.time_code(rt['wt_galgal_counts'], unit='min') #.TE. replace start time with runtime in minutes


    # random random
    autocorr=1
    print('Calculating RR_counts...')
    rt['wt_randrand_counts'] = hf.time_code('start') #.TS. get code start time
    RR_counts = DDtheta_mocks(autocorr, nthreads, tbins, rand_RA, rand_DEC) #.TC.
    rt['wt_randrand_counts'] = hf.time_code(rt['wt_randrand_counts'], unit='min') #.TE. replace start time with runtime in minutes


    # gal random
    autocorr=0
    print('Calculating DR_counts...')
    rt['wt_galrand_counts'] = hf.time_code('start') #.TS. get code start time
    DR_counts = DDtheta_mocks(autocorr, nthreads, tbins, RA, DEC, RA2=rand_RA, DEC2=rand_DEC) #.TC.
    rt['wt_galrand_counts'] = hf.time_code(rt['wt_galrand_counts'], unit='min') #.TE. replace start time with runtime in minutes
    #---

    # calc w(theta)
    print('Converting counts to correlation...')
    rt['wt_counts_to_cf'] = hf.time_code('start') #.TS. get code start time
    wtheta = convert_3d_counts_to_cf(N, N, rand_N, rand_N, DD_counts, DR_counts, DR_counts, RR_counts) #.TC.
    rt['wt_counts_to_cf'] = hf.time_code(rt['wt_counts_to_cf'], unit='min') #.TE. replace start time with runtime in minutes

    bcens = np.around((tbins[:-1]+tbins[1:])/2, decimals=5)
    return [bcens, wtheta, rt]