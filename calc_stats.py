import numpy as np
import pandas as pd
from pathlib import Path

from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks
from Corrfunc.utils import convert_3d_counts_to_cf
from Corrfunc.theory.wp import wp
from Corrfunc.theory.xi import xi


def calc_stats(param_dict, boxes):
    """
    Args:
        boxes (dict): keys {gals_PS', 'gals_RDZ', 'rands_PS', 'rands_RDZ'}
    """
    p = param_dict

    ### Calc stats for each redshift bin
    # wtheta
    if 'wtheta' in p['stats']:
        # precision for get_group('zbin')
        zbins = list(boxes['gals_RDZ'].groupby('zbin').groups.keys())
        zbins_dict = dict(zip(np.round(zbins, 2), zbins))

        for zzz in p['zbin_cens']: # calculate for each zbin
            theta_binedges, rp_binedges = get_theta_rp_from_tratio_bins(p, zzz)
            theta_bincens = calc_bincens(theta_binedges, decimals=5)

            zz = zbins_dict[np.round(zzz,2)]
            gals_RDZ_z = boxes['gals_RDZ'].groupby('zbin').get_group(zz)
            rands_RDZ_z = boxes['rands_RDZ'].groupby('zbin').get_group(zz)

            wtheta = calc_wtheta(p, gals_RDZ_z, rands_RDZ_z, theta_binedges)
            write_stat_to_file( p, 'wtheta', wtheta, theta_bincens, zzz, p['zbin_width'], \
                                len(gals_RDZ_z), len(rands_RDZ_z) )

    # theory stats
    gals_PS = boxes['gals_PS']
    numgals = len(gals_PS.index)
    zmid = np.round((gals_PS.Redshift.max()+gals_PS.Redshift.min())/2, 2)
    zwidth = np.round(gals_PS.Redshift.max() - gals_PS.Redshift.min(), 2)
    rbcens = calc_bincens(p['rbin_edges'], decimals=5)

    if 'xi' in p['stats']:
        xi = calc_xi(p, gals_PS)
        write_stat_to_file( p, 'xi', xi, rbcens, zmid, zwidth, numgals, np.nan )

    if 'wp' in p['stats']:
        wp = calc_wp(p, gals_PS)
        write_stat_to_file( p, 'wp', wp, rbcens, zmid, zwidth, numgals, np.nan )
    ###

    return None



def get_theta_rp_from_tratio_bins(param_dict, redshift, invert=False):
    """ Converts tratio_binedges from units of theta_BAO(zbin) to degrees
        at given redshift using law of cosines.

        invert == True: 2nd arg (tratio_binedges) should be t_binedges
                        Returns tratio_binedges
    """

    p = param_dict

    d_BAO = 105 # h^-1 Mpc
    rz = (p['cosmo'].comoving_distance(redshift).value) * p['cosmo'].h # Mpc/h
    theta_BAO = np.degrees(np.arccos(1 - d_BAO**2/2/rz**2)) # law of cosines

    if not invert:
        t_binedges = theta_BAO * p['theta_scaled_binedges'] # degrees
        rp_binedges = rz * np.sqrt(2*(1-np.cos(np.radians(t_binedges)))) # projected, comoving dist
        return t_binedges, rp_binedges
    else:
        trat_binedges = p['theta_scaled_binedges'] / theta_BAO
        return trat_binedges

    return None


def calc_bincens(bin_edges, decimals=5):
    """ Returns centers of bins defined by bin_edges.
        Array length = len(bin_edges) - 1.
    """
    return np.around((bin_edges[:-1]+bin_edges[1:])/2, decimals=decimals)


def calc_wtheta(param_dict, gals_RDZ, rands_RDZ, tbins):
    """
        gals_RDZ = galaxies DataFrame including (at least) columns 'RA' and 'DEC'
        rands_RDZ = randoms DataFrame including (at least) columns 'RA' and 'DEC'
        tbins (array): theta bin edges [degrees]

    Returns wtheta in each bin
    """
    p = param_dict
    nthreads = p['nthreads']

    RA, DEC, N = np.asarray(gals_RDZ.RA), np.asarray(gals_RDZ.DEC), len(gals_RDZ)
    ran_RA, ran_DEC, ran_N = np.asarray(rands_RDZ.RA), np.asarray(rands_RDZ.DEC), len(rands_RDZ)

    #--- Calc pair counts---

    # DDtheta_mocks expects ra in [0.0, 360.0] and dec in [-90.0, 90.0] degrees
    # DDtheta_mocks returns numpy structured array containing
    #   [thetamin, thetamax, thetaavg, npairs, weightavg] for each angular bin

    # gal gal
    autocorr=1
    print('Calculating DD_cnts...')
    DD_cnts = DDtheta_mocks(autocorr, nthreads, tbins, RA, DEC)


    # random random
    autocorr=1
    print('Calculating RR_cnts...')
    RR_cnts = DDtheta_mocks(autocorr, nthreads, tbins, ran_RA, ran_DEC)


    # gal random
    autocorr=0
    print('Calculating DR_cnts...')
    DR_cnts = DDtheta_mocks(autocorr, nthreads, tbins, RA, DEC, RA2=ran_RA, DEC2=ran_DEC)
    #---

    # calc w(theta)
    print('Converting counts to correlation...')
    wtheta = convert_3d_counts_to_cf(N, N, ran_N, ran_N, DD_cnts, DR_cnts, DR_cnts, RR_cnts)

    return wtheta


def calc_xi(param_dict, gals_PS):
    p = param_dict
    nthreads = p['nthreads']
    boxsize = p['mock_Lbox']
    bins = p['rbin_edges']

    X,Y,Z = gals_PS.x_theory, gals_PS.y_theory, gals_PS.z_theory
    results = xi(boxsize, nthreads, bins, X, Y, Z)

    return results['xi']


def calc_wp(param_dict, gals_PS):
    p = param_dict
    nthreads = p['nthreads']
    boxsize = p['mock_Lbox']
    bins = p['rbin_edges']
    pimax = p['pimax']

    X,Y,Z = gals_PS.x_theory, gals_PS.y_theory, gals_PS.z_theory
    results = wp(boxsize, pimax, nthreads, bins, X, Y, Z)

    return results['wp']


def write_stat_to_file(param_dict, statname, statdat, bincens, zbin, zwidth, Ngals, Nrands):
    """ statname (string): name of statistic
        statdat (array): stat calculated in each bin. must be same length as bincens
        bincens (array): center of each bin for which the stat was calculated.
        zbin = center of the redshift bin for which the stat was calculated.

        Checks that existing file (if any) has same number of columns
    """
    p = param_dict
    numbins = len(bincens)
    colwidth = 25
    # if you change this you must update several things in the rest of this function:
    extra_cols = ['mock', 'zbin', 'zwidth', 'Nstack', 'Ngals', 'Nrands', 'statname']
    lc = len(extra_cols) + 2*numbins # num cols to write to file

    # Check if file exists and has same number of columns as expected
    fpath = Path(p['statfout'])
    if fpath.is_file():
        try:
            df = pd.read_csv(fpath, comment='#', delim_whitespace=True, nrows=10)
            # get the structure of the current file
            lfc = len(df.columns) # num cols in existing file
            assert lfc==lc
        except: # if df can't be read or lengths don't match, move existing
            mv_fout = hf.file_ow(p['statfout'])
            print('*** Stat file format incompatible.')
            print('\tMoved existing file to {} so it is not overwritten. ***'.format(mv_fout))

    # If statfout has been moved (above) or never existed, create new file.
    if not fpath.is_file():
        hdr =  ('bin_[#] columns contain bin centers (theta in degrees, r in Mpc/h). '
                'stat_[#] contain the stat for each bin, others are extra info.\n')
        bcols = ['bin_{}'.format(i) for i in range(numbins)]
        scols = ['stat_{}'.format(i) for i in range(numbins)]
        new_cols = ''.join(i.rjust(colwidth) for i in extra_cols+bcols+scols)
        print('# {}\n{}'.format(hdr, new_cols), file=open(p['statfout'], 'w'))
        # np.savetxt(self.statfout, new_cols, header=hdr) # write header

    # Now fout exists, so append the data.
    print("Appending stat '{}' to {}".format(statname, p['statfout']))
    # mlw = 50*lc
    dat_xcols = ['{}'.format(p['mock_num']), '{}'.format(zbin), '{}'.format(zwidth),
                 '{}'.format(p['mock_Nstack']), '{:.5e}'.format(Ngals), '{:.5e}'.format(Nrands),
                 statname ]
    dat_bcols = ['{:.15f}'.format(bincens[b]) for b in range(numbins)]
    dat_scols = ['{:.15f}'.format(statdat[s]) for s in range(numbins)]
    str_cols = ''.join(i.rjust(colwidth) for i in dat_xcols+dat_bcols+dat_scols)
    # str_cols = np.array2string(np.append(dat_xcols, np.append(bincens, statdat)), \
    #         formatter={'float_kind':lambda x: "%25.15e" % x, 'str':lambda x: "%25s" % x},
     # max_line_width=mlw)[1:-1]
    print(str_cols, file=open(p['statfout'], 'a'))

    return None
