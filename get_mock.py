""" Generates a mock galaxy catalog by loading a DM simulation, identifying halos
    if necessary, and populating with galaxies using an HOD.
    Also returns a box of randoms at the same coordinates (optional).
"""

import halotools.sim_manager.sim_defaults as sim_defaults
from halotools.sim_manager import CachedHaloCatalog
from halotools.empirical_models import PrebuiltHodModelFactory


def get_mock(param_dict, randoms=True):
    """ Returns a DF of galaxies populated on a DM sim and a box of randoms (optional).

    Args:

        param_dict  (dict): Must have keys: 'sim_name',

        randoms     (bool): If True, returns a box of randoms with the same coordinate bounds
                            as the galaxy box.

    Returns:

        gals_xyz   (dataframe): Mock box with rows of galaxy coordinates

        rands_xyz  (dataframe): Mock box with rows of randoms coordinates
    """

    p = param_dict
    load_sim = { # dict of functions to load a DM simulation and return halo coords
                'multidark': load_multidark,
                'outerrim': load_outerrim
                }

    halo_xyz = load_sim[p['sim_name']](p)
    gals_xyz = popGals_usingHOD(halo_xyz, p)

    return None


def load_multidark(param_dict):
    """ Loads cached Multidark particles using Halotools
    """
    p = param_dict

    sim_defaults.default_simname = p['sim_name']
    sim_defaults.default_redshift = p['sim_redshift'] # multidark actual z = 0.466

    # get Multidark + ROCKSTAR halo catalog
    halo_xyz = CachedHaloCatalog(simname=p['sim_name'], halo_finder=p['sim_halofinder'], \
                version_name = 'halotools_v0p4', redshift=p['sim_redshift'])

    return halo_xyz


def load_outerrim(param_dict):

    return halo_xyz


def popGals_usingHOD(halo_xyz, param_dict):

    HODmodel = PrebuiltHodModelFactory(param_dict['HOD_model'], redshift=redshift)
    HODmodel.param_dict['logMmin'] = HODparams[0] # Minimum mass required for a halo to host a central galaxy.
    HODmodel.param_dict['sigma_logM'] = HODparams[1] # Rate of transition from ⟨Ncen⟩=0⇒⟨Ncen=1⟩
    HODmodel.param_dict['logM0'] = HODparams[2] # Low-mass cutoff in ⟨Nsat⟩
    HODmodel.param_dict['logM1'] = HODparams[3] # Characteristic halo mass where ⟨Nsat⟩ begins to assume a power law form.
    HODmodel.param_dict['alpha'] = HODparams[4] # Power law slope of the relation between halo mass and ⟨Nsat⟩.
    HODmodel.param_dict
