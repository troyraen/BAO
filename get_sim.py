""" Returns a mock galaxy catalog (phase space coords) by loading a DM simulation,
    identifying halos if necessary, and populating with galaxies using an HOD.

    Use Arg: param_dict = main.load_param_dict()
"""

import os as os
import pandas as pd

import halotools.sim_manager.sim_defaults as sim_defaults
from halotools.sim_manager import CachedHaloCatalog
from halotools.sim_manager import UserSuppliedHaloCatalog
from halotools.empirical_models import PrebuiltHodModelFactory

import generic_io

def get_sim_galaxies(param_dict):
    """ Returns a DF of galaxies populated on a DM sim.

    Args:
        param_dict  (dict): Must have keys: 'sim_name',

    Returns:
        gals_PS   (dataframe): Mock box with rows of galaxy phase space coordinates
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

    p = param_dict

    # Load halo data from Outer Rim files
    halo_files, read_cols, rename_cols = load_outerrim_data_setup(p)
    coldat = []
    for hf in halo_files:
        gio = generic_io.Generic_IO(hf, None)
        coldat.append(gio.read_columns(read_cols))
    halo_df = pd.concat(coldat, axis=0).rename(columns = rename_cols)

    # Generate Halotools halo catalog
    metadata, halodata = load_outerrim_halotools_setup(p, halo_df)
    halocat = UserSuppliedHaloCatalog(**metadata, **halodata)

    return halocat


def load_outerrim_data_setup(param_dict):
    """
    Returns:

        halo_files  (list): paths to all data files at redshift=param_dict['sim_redshift']

        read_cols   (list): list of file columns to read in

        rename_cols (dict): columns to rename for Halotools input
    """

    p = param_dict

    # Directory structure
    halo_dir = "/home/tjr63/Osiris/BAO_simdata/OuterRim/M000/L4225/HACC000/analysis/" + \
                "Halos/b0168/FOFp/"
    halo_files_ztail = { # key = redshift, value = list of files
                    0.502242:
                        ['STEP331/m000-331.fofproperties'] + \
                        ['STEP331/m000-331.fofproperties#{}'.format(n) for n in range(87,105)],
                    0.539051:
                        ['STEP323/m000-323.fofproperties'] + \
                        ['STEP323/m000-323.fofproperties#{}'.format(n) for n in range(87,105)]
                  }
    halo_files = [ os.path.join(halo_dir,hf) for hf in halo_files_ztail[p['sim_redshift']] ]

    # Columns to read and rename
    read_cols = [ 'fof_halo_count', 'fof_halo_tag', 'fof_halo_mass', \
                  'fof_halo_center_x', 'fof_halo_center_y', 'fof_halo_center_z' ]

    rename_cols = { 'fof_halo_tag': 'halo_id',
                    'fof_halo_mass': 'halo_mass',
                    'fof_halo_count': 'halo_count',
                    'fof_halo_center_x': 'halo_x',
                    'fof_halo_center_y': 'halo_y',
                    'fof_halo_center_z': 'halo_z' }

    return halo_files, read_cols, rename_cols


def load_outerrim_halotools_setup(param_dict, halo_df):
    """ https://halotools.readthedocs.io/en/latest/api/halotools.sim_manager
            .UserSuppliedHaloCatalog.html#halotools.sim_manager.UserSuppliedHaloCatalog
    """

    p = param_dict

    metadata = {'Lbox': p['mock_Lbox'],
                'particle_mass': p['sim_particle_mass'],
                'redshift': p['z4push']
                }

    halodata = halo_df.to_dict(orient='series')

    return metadata, halodata


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
