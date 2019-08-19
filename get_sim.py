""" Returns a mock galaxy catalog (phase space coords) by loading a DM simulation,
    identifying halos if necessary, and populating with galaxies using an HOD.
    # Also returns a box of randoms at the same coordinates (optional).

    Use Arg: param_dict = main.load_param_dict()
"""

import halotools.sim_manager.sim_defaults as sim_defaults
from halotools.sim_manager import CachedHaloCatalog
from halotools.empirical_models import PrebuiltHodModelFactory


def get_sim_galaxies(param_dict, randoms=True):
    """ Returns a DF of galaxies populated on a DM sim and a box of randoms (optional).

    Args:

        param_dict  (dict): Must have keys: 'sim_name',

        randoms     (bool): If True, returns a box of randoms with the same coordinate bounds
                            as the galaxy box.

    Returns:

        gals_PS   (dataframe): Mock box with rows of galaxy phase space coordinates

        rands_PS  (dataframe): Mock box with rows of randoms phase space coordinates
    """

    p = param_dict
    load_sim = { # dict of fncs to load a DM sim and return a Halotools halo_catalog object
                'multidark': load_multidark,
                'outerrim': load_outerrim
                }

    halocat = load_sim[p['sim_name']](p)
    gals_PS = popHalos_usingHOD(halocat, p) # phase space coords of galaxies

    # check and set param_dict parameters
    assert halocat.Lbox[0] == p['sim_Lbox'], 'sim_Lbox != cat_Lbox'
    assert halocat.redshift == p['sim_redshift'], 'sim_redshift != cat_redshift'

    return gals_PS


def load_multidark(param_dict):
    """ Loads cached Multidark particles using Halotools.
        Returns a Halotools halo_catalog object.
    """
    p = param_dict

    sim_defaults.default_simname = p['sim_name']
    sim_defaults.default_redshift = p['sim_redshift'] # multidark actual z = 0.466

    # get Multidark + ROCKSTAR halo catalog
    halocat = CachedHaloCatalog(simname=p['sim_name'], halo_finder=p['sim_halofinder'], \
                version_name = 'halotools_v0p4', redshift=p['sim_redshift'])

    return halocat


def load_outerrim(param_dict):

    return halocat


def popHalos_usingHOD(halocat, param_dict):
    """ Populates a Halotools halo catalog with galaxies using an HOD.
        Returns a DataFrame of galaxy phase space coordinates.
    """

    p = param_dict

    # Populate the catalog
    HODmodel = PrebuiltHodModelFactory(p['HOD_model'], redshift=p['sim_redshift'])
    for key, val in p['HOD_params'].items():
        HODmodel.param_dict[key] = val

    HODmodel.populate_mock(halocat) # use this command to create and populate the mock
    # HODmodel.mock.populate() # use this command to REpopulate

    # Convert to DF
    gals = HODmodel.mock.galaxy_table.to_pandas()

    return gals[['x','y','z', 'vx','vy','vz']]


# def get_randoms(param_dict):
#
#     p = param_dict
#
#     # create random points in box with side length self.Lbox, centered around origin
#     ran_coords = np.random.random((self.Nrands,3))*self.Lbox - self.Lbox/2
#     ran_vels = np.zeros((self.Nrands,3))
#     # Set self.RandomsRDZ DF
#     self.RandomsPS = pd.DataFrame(np.hstack([ran_coords,ran_vels]), columns=['x','y','z', 'vx','vy','vz'])
