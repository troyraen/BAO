import numpy as np
import pandas as pd
from astropy import cosmology

import get_sim as gs
import transform_sim as ts
import plots as plots


#*** param_dict keys ***#
pdkeys_calc = [     # parameters that need to be calculated
                    'cosmo', # astropy cosmology object
                    'mock_Lbox',
                    'theta_scaled_binedges',
                    'zbin_cens', # set in transform_sim.transform()
                    'zbin_edges'  # set in transform_sim.transform()
                    ]
pdkeys_noncalc = [  # parameters set as is
                    'cosmo_H0',
                    'cosmo_Om0',
                    'galplots',
                    'HOD_model',
                    'HOD_params',
                    'Nrands',
                    'mock_Nstack', # num sim boxes to stack for mock.  1 or even integer.
                    'pimax',
                    'sim_halofinder',
                    'sim_Lbox',
                    'sim_name',
                    'sim_redshift',
                    'stats',
                    'statfout',
                    'theta_scaled_min',
                    'theta_scaled_max',
                    'theta_scaled_Nbins',
                    'z4push',
                    'zbin_width',
                    ]
pdkeys = pdkeys_noncalc + pdkeys_calc


# fs*** param_dict DEFAULTS ***#

# simulation info
sim_Lbox = 1000.0 # {'multidark':1000.0, 'outerrim':}
sim_name = 'outerrim' # which DM simulation to use
sim_redshift = 0.466 # {'multidark':0.466, 'outerrim':}
sim_halofinder = 'rockstar' # outerrim loads halos directly (found using FoF)

# cosmology, HOD info
cosmo_H0 = 70.0
cosmo_Om0 = 0.3
HOD_model = 'zheng07'
# https://halotools.readthedocs.io/en/latest/quickstart_and_tutorials/
# tutorials/model_building/preloaded_models/zheng07_composite_model.html
HOD_params = {  'logMmin': 12.89376701, # Minimum mass req for halo to host central galaxy.
                'sigma_logM': 0.23939566, # Rate of transition from ⟨Ncen⟩=0⇒⟨Ncen=1⟩
                'logM0': 12.26827089, # Low-mass cutoff in ⟨Nsat⟩
                'logM1': 14.03372441, # Halo mass where ⟨Nsat⟩ begins to assume power law form.
                'alpha': 1.32828278 # Power law slope of halo mass, ⟨Nsat⟩ relation.
            }

# mock info
z4push = 'sim' # either a float or 'sim'. start box coords at this redshift distance
zbin_width = 0.3 # redshift bin width
# theta bins
# theta_scaled_: theta bins in units of the BAO scale
theta_scaled_min = 0.05 # degrees/theta_bao. min separation to calc wtheta
theta_scaled_max = 2. # degrees/theta_bao. max separation to calc wtheta
theta_scaled_Nbins = 50 # number of theta bins for wtheta calculation

# misc
galplots = False # whether to plot galaxies while obtaining/transforming
mock_Nstack = 1
Nrands = 5
pimax = 300
stats = ['wtheta']#, 'xi', 'wp'] # which stats to calculate
statfout = 'data/stats.dat' # file path to write stats

# fe*** END param_dict DEFAULTS ***#


def proc_mockbox(param_dict={}):
    """ Loads DM sim particles or halos.
        Populates sim with galaxies.
        Transforms coordinates to ra/dec/redshift with box at redshift=z4push.
        Calculates correlation statistics and writes them to file at statfout.

    Args:

        param_dict  (dict): Parameters for the run.
                            Keys can be anything in pdkeys_noncalc (defined above).
                            Default values (defined above) are used for missing keys.
    """
    # Load parameter dictionary for the run
    pdict = load_param_dict(param_dict)
    # print(pdict)

    # Load galaxy box from DM sim
    gals_PS = gs.get_sim_galaxies(pdict, randoms=False)
    if pdict['galplots']:
        plots.plot_galaxies(gals_PS, title="Sim Galaxies")
    # print(gals_PS.sample(2))

    # Transform coordinates
    gals_PS, gals_RDZ, pdict = ts.transform(pdict, gals_PS)
    if pdict['galplots']:
        plots.plot_galaxies(gals_PS, title="Sim Galaxies Transformed")
        plots.plot_galaxies(gals_PS, title="Sim Galaxies Transformed", coords='rz')
    # print(pdict)
    # print(gals_PS.sample(2))
    # print(gals_RDZ.sample(2))


    # Get randoms
    rands_PS, rands_RDZ = get_randoms(pdict)
    if pdict['galplots']:
        plots.plot_galaxies(rands_PS, title="Randoms")
        plots.plot_galaxies(rands_PS, title="Randoms", coords='rz')
    # print(rands_PS.sample(2))
    # print(rands_RDZ.sample(2))

    return None


def load_param_dict(param_dict={}):
    """ Loads parameter dictionary for the run.

    Args:

        param_dict  (dict): Parameters for the run.
                            Keys can be anything in pdkeys_noncalc (defined above).
                            Default values (defined above) are used for missing keys.
    """
    pdict = param_dict.copy() # set default and calculated values below.

    # Set non-calculated defaults for keys not already present
    for p in (set(pdkeys_noncalc) - set(param_dict.keys())):
        pdict[p] = globals()[p]

    # Set calculated parameters
    pdict['cosmo'] = cosmology.FlatLambdaCDM(H0=pdict['cosmo_H0'],
                                             Om0=pdict['cosmo_Om0'])

    pdict['mock_Lbox'] = pdict['sim_Lbox']* pdict['mock_Nstack']

    pdict['theta_scaled_binedges'] = np.logspace(np.log10(pdict['theta_scaled_min']),
                                                 np.log10(pdict['theta_scaled_max']),
                                                 pdict['theta_scaled_Nbins']+1),

    return pdict


def get_randoms(param_dict):

    p = param_dict

    # create random points in box, centered around origin
    Lbox, Nrands = p['mock_Lbox'], p['Nrands']
    ran_coords = np.random.random((Nrands,3))*Lbox - Lbox/2
    ran_vels = np.zeros((Nrands,3))
    rands_PS = pd.DataFrame(np.hstack([ran_coords,ran_vels]),
                                        columns=['x','y','z', 'vx','vy','vz'])

    rands_PS = ts.shift_face_to_z4push(p, rands_PS)
    rands_RDZ, _,_ = ts.get_ra_dec_z(p, rands_PS)

    return rands_PS, rands_RDZ


                # mb = MBox( Nstack=Nstack, z4push=z4push, zbin_width=zw, \
                #             tratio_binedges=tratio_binedges, rbin_edges=rbin_edges, \
                #             pimax=pimax, Nrands=Nrands, galplots=galplots, \
                #             statfout=statfout
                #         )
                # mb.getmock(simname=simname)
                # mb.calc_stats(stats=stats, nthreads=24)
