# fs*** IMPORTS ***#
import numpy as np
import pandas as pd
from pathlib import Path
import datetime
from astropy import cosmology

import get_sim as gs
import transform_sim as ts
import plots as plots
import calc_stats as cs
# fe*** IMPORTS ***#

#*** param_dict keys ***#
pdkeys_calc = [     # parameters that need to be calculated
                    'cosmo', # astropy cosmology object
                    'mock_Lbox',
                    'mock_num', # ID for this mock
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
                    'Nrands', # int or 'sim'. If 'sim', Nrands set in transform_sim.transform()
                    'nthreads', # num threads to use in corr_func calculations.
                    'mock_Nstack', # num sim boxes to stack for mock.  1 or even integer.
                    'pimax',
                    'rbin_edges',
                    'sim_halofinder',
                    'sim_Lbox', # h^-1 [Mpc]
                    'sim_name',
                    'sim_redshift',
                    'stats',
                    'statfout',
                    'theta_scaled_min',
                    'theta_scaled_max',
                    'theta_scaled_Nbins',
                    'z4push', # float or 'sim'. if 'sim', z4push set in transform_sim.transform()
                    'zbin_width',
                    ]
pdkeys = pdkeys_noncalc + pdkeys_calc


# fs*** param_dict DEFAULTS ***#

# Simulation info
sim_Lbox = 1000.0 # {'multidark':1000.0, 'outerrim':}
sim_name = 'multidark' # which DM simulation to use
sim_redshift = 0.466 # {'multidark':0.466, 'outerrim':}
sim_halofinder = 'rockstar' # outerrim loads halos directly (found using FoF)

# Cosmology, HOD info
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

# Mock info
mock_Nstack = 1
Nrands = 'sim'
z4push = 'sim'

# Stats info
nthreads = 24
pimax = 300
rbin_edges = np.logspace(np.log10(25.0), np.log10(150.0), 51)
stats = ['wtheta']#, 'xi', 'wp'] # which stats to calculate
statfout = Path(__file__).resolve().parent / 'data/stats.dat' # file to write stats
# theta bins. theta_scaled_ => theta bins in units of the BAO scale
theta_scaled_min = 0.05 # degrees/theta_bao. min separation to calc wtheta
theta_scaled_max = 2. # degrees/theta_bao. max separation to calc wtheta
theta_scaled_Nbins = 50 # number of theta bins for wtheta calculation
zbin_width = 0.3 # wtheta redshift bins

# Misc
galplots = False # whether to plot galaxies while obtaining/transforming

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
    Returns:
        p['statfout'] (path): path to stats output.
    """
    # Load parameter dictionary for the run
    p = load_param_dict(param_dict)
    # print(p)

    # Load galaxy box from DM sim
    gals_PS = gs.get_sim_galaxies(p)
    if p['galplots']: plots.plot_galaxies(gals_PS, title="Sim Galaxies")
    # print(gals_PS.info())

    # Transform coordinates
    gals_PS, gals_RDZ, p = ts.transform(p, gals_PS)
    if p['galplots']:
        plots.plot_galaxies(gals_PS, title="Sim Galaxies Transformed")
        plots.plot_galaxies(gals_PS, title="Sim Galaxies Transformed", coords='rz')
    # print(p)
    # print(gals_PS.sample(2))
    # print(gals_RDZ.sample(2))

    # Get randoms
    rands_PS, rands_RDZ = get_randoms(p)
    # print(rands_PS.info())
    if p['galplots']:
        plots.plot_galaxies(rands_PS, title="Randoms")
        plots.plot_galaxies(rands_PS, title="Randoms", coords='rz')
    # print(rands_PS.sample(2))
    # print(rands_RDZ.sample(2))


    # Calc stats
    boxes = { 'gals_PS':gals_PS, 'gals_RDZ':gals_RDZ,
              'rands_PS':rands_PS, 'rands_RDZ':rands_RDZ }
    cs.calc_stats(p, boxes)

    return p['statfout']


def load_param_dict(param_dict={}):
    """ Loads parameter dictionary for the run.

    Args:

        param_dict  (dict): Parameters for the run.
                            Keys can be anything in pdkeys_noncalc (defined above).
                            Default values (defined above) are used for missing keys.
    """
    p = param_dict.copy() # set default and calculated values below.

    # Set non-calculated defaults for keys not already present
    for k in (set(pdkeys_noncalc) - set(param_dict.keys())):
        p[k] = globals()[k]

    # Set calculated parameters
    p['cosmo'] = cosmology.FlatLambdaCDM(H0=p['cosmo_H0'],
                                             Om0=p['cosmo_Om0'])

    p['mock_Lbox'] = p['sim_Lbox']* p['mock_Nstack']

    p['mock_num'] = get_mock_num()

    p['theta_scaled_binedges'] = np.logspace(np.log10(p['theta_scaled_min']),
                                                 np.log10(p['theta_scaled_max']),
                                                 p['theta_scaled_Nbins']+1)

    return p


def get_mock_num():
    dtm = datetime.datetime.now() # get date and time to use as mock number
    mocknum = float(dtm.strftime("%m%d%y.%H%M"))
    return mocknum


def get_randoms(param_dict):

    p = param_dict

    # create random points in box, centered around origin
    Lbox, Nrands = p['mock_Lbox'], int(p['Nrands'])
    ran_coords = np.random.random((Nrands,3))*Lbox - Lbox/2
    ran_vels = np.zeros((Nrands,3))
    rands_PS = pd.DataFrame(np.hstack([ran_coords,ran_vels]),
                                        columns=['x','y','z', 'vx','vy','vz'])

    rands_PS, _ = ts.shift_face_to_z4push(p, rands_PS)
    rands_PS, rands_RDZ, _,_ = ts.get_ra_dec_z(p, rands_PS)

    return rands_PS, rands_RDZ
