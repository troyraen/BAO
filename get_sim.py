""" Returns a mock galaxy catalog (phase space coords) by loading a DM simulation,
    identifying halos if necessary, and populating with galaxies using an HOD.

    Use Arg: param_dict = main.load_param_dict()
"""

import os as os
from warnings import warn as _warn
import random as rand
import pandas as pd
import numpy as np

import halotools.sim_manager.sim_defaults as sim_defaults
from halotools.sim_manager import CachedHaloCatalog
from halotools.sim_manager import UserSuppliedHaloCatalog
from halotools.empirical_models import PrebuiltHodModelFactory

import genericio as gio
# import generic_io

def get_sim_galaxies(param_dict):
    """ Returns a DF of galaxies populated on a DM sim.

    Args:
        param_dict  (dict): as returned by main.load_param_dict()

    Returns:
        gals_PS   (dataframe): Mock box with rows of galaxy phase space coordinates
    """

    p = param_dict
    load_sim = { # dict of fncs to load a DM sim and return a Halotools halo_catalog object
                'multidark': load_multidark,
                'outerrim': load_outerrim,
                'outerrim_gioforked': load_outerrim_gioforked
                }

    halocat = load_sim[p['sim_name']](p)
    gals_PS = popHalos_usingHOD(halocat, p) # phase space coords of galaxies

    # check and set param_dict parameters
    # these are probably unnecessary
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
    halo_metafile, read_cols, name_df_cols = load_outerrim_data_setup(p)
    print('Loading Outer Rim halos from file.')
    data = gio.read(halo_metafile, read_cols)  # nparray [len(read_cols), num halos]

    # Generate Halotools halo catalog
    metadata, halodata = load_outerrim_halotools_setup(p, data, name_df_cols)
    del data
    halocat = UserSuppliedHaloCatalog(**metadata, **halodata)

    return halocat


def load_outerrim_data_setup(param_dict):
    """
    Returns:

        halo_metafile (str): path to non-hashed Outer Rim halos file at
                             redshift = param_dict['sim_redshift']

        read_cols    (list): list of file columns to read in

        name_df_cols (list): renamed columns for Halotools input
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
    halo_metafile = halo_files[0]

    # Columns to read and rename
    read_cols = [ 'fof_halo_tag', 'fof_halo_mass', \
                  'fof_halo_center_x', 'fof_halo_center_y', 'fof_halo_center_z',
                  'fof_halo_mean_vx', 'fof_halo_mean_vy', 'fof_halo_mean_vz' ]

    name_df_cols = [ 'halo_id', 'halo_mvir', 'halo_x', 'halo_y', 'halo_z',
                     'halo_vx', 'halo_vy', 'halo_vz' ]

    return halo_metafile, read_cols, name_df_cols


def load_outerrim_halotools_setup(param_dict, data, name_df_cols):
    """ https://halotools.readthedocs.io/en/latest/api/halotools.sim_manager
            .UserSuppliedHaloCatalog.html#halotools.sim_manager.UserSuppliedHaloCatalog

    Args:
        param_dict as returned by main.load_param_dict()

        data      (nparray): as returned by gio.read()

        name_df_cols (list): as returned by load_outerrim_data_setup()
    """

    p = param_dict

    metadata = {'Lbox': p['sim_Lbox'],
                'particle_mass': p['sim_particle_mass'],
                'redshift': p['sim_redshift']
                }

    # create df with halos having x < 500 Mpc/h
    # _warn('\nThrowing out all halos with x > 500 Mpc/h\n')
    # halo_df = pd.DataFrame(data.T[data.T[:,2]<500,:], columns=name_df_cols)
    halofrac = p['keep_halo_frac']
    if halofrac != 1:
        _warn("\nDownsampling halos\n")
        l = data.shape[1]
        keep_idx = rand.sample(range(l),int(l*halofrac))
        halo_df = pd.DataFrame(data.T[keep_idx,:], columns=name_df_cols)
    else:
        halo_df = pd.DataFrame(data.T, columns=name_df_cols)

    # add columns
    halo_df['halo_rvir'] = load_outerrim_halotools_setup_calc_rvir(p, halo_df)
                           # to calc concentration
    halo_df['halo_upid'] = -1 # indicate all are host halos
    halo_df['halo_hostid'] = halo_df['halo_id']

    # get dict of arrays for halotools
    halodata = halo_df.to_dict(orient='series')
    for k, val in halodata.items():
        halodata[k] = val.to_numpy()

    return metadata, halodata


def load_outerrim_halotools_setup_calc_rvir(param_dict, halo_df):
    """ Calculates and returns an estimated r_vir for each halo.
        Assumes fof_halo_mass = M_vir
            (column name already set from load_outerrim_data_setup()).
    """
    p = param_dict

    # estimate concentration using concentration-mass relation (Dutton+14, eqs 7, 12, 13)
    z = p['sim_redshift']
    conc_A = 0.537 + (1.025-0.537) * np.exp(-0.718 * z**1.08)
    conc_B = -0.097 + 0.024 * z
    conc = 10**conc_A * (halo_df.halo_mvir/10**12)**conc_B

    # estimate overdensity using linking length and concentration (More+11 eq 13)
    cp1 = conc+1
    mu = np.log(cp1) - conc/cp1
    psi = conc**2/mu/cp1**2
    Delta_vir = (p['sim_FoF_b']/0.2)**(-3) * 244.86/psi - 1

    # calculate r_vir
    Mpc = 3.086e22 # meters
    G = 6.674e-11 # m^3/kg/s^2
    Msun = 1.989e30 # kg
    rho_mean = 3e10*p['cosmo_Om0']*Mpc/(8*np.pi*G) # h^2 kg Mpc^-3
    rvir = ((3*halo_df.halo_mvir*Msun)/(4*np.pi* Delta_vir* rho_mean))**(1./3.) # Mpc/h

    return rvir # todo: check this better


def popHalos_usingHOD(halocat, param_dict):
    """ Populates a Halotools halo catalog with galaxies using an HOD.
        Returns a DataFrame of galaxy phase space coordinates.
    """

    p = param_dict
    simz = p['sim_redshift']

    # Load the HOD
    kwargs = {'redshift': simz}
    if p['sim_name'] == 'outerrim': # need to calc concentration
        kwargs['conc_mass_model'] = 'dutton_maccio14'
    HODmodel = PrebuiltHodModelFactory(p['HOD_model'], **kwargs)
    # add HOD params
    for key, hoddic in p['HOD_params'].items():
        if simz>key[0] and simz<=key[1]: # get the right params for sim redshift
            for key, val in hoddic.items():
                HODmodel.param_dict[key] = val
        break

    # Populate the catalog
    HODmodel.populate_mock(halocat) # use this command to create and populate the mock
    # HODmodel.mock.populate() # use this command to REpopulate

    # Convert to DF
    gals = HODmodel.mock.galaxy_table.to_pandas()

    return gals[['x','y','z', 'vx','vy','vz']]


def load_outerrim_gioforked(param_dict):

    p = param_dict

    # Load halo data from Outer Rim files
    halo_files, read_cols, rename_cols = load_outerrim_gioforked_data_setup(p)
    coldat = []
    for hf in halo_files:
        gio = generic_io.Generic_IO(hf, None)
        coldat.append(gio.read_columns(read_cols))
    halo_df = pd.concat(coldat, axis=0).rename(columns = rename_cols)

    # Add column indicating host halo
    halo_df['halo_upid'] = -1

    # Generate Halotools halo catalog
    metadata, halodata = load_outerrim_halotools_setup(p, halo_df)
    halocat = UserSuppliedHaloCatalog(**metadata, **halodata)

    return halocat


def load_outerrim_gioforked_data_setup(param_dict):
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
