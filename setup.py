import numpy as np
import pandas as pd

from astropy import cosmology
# from astropy.constants import c  # the speed of light

import halotools.sim_manager.sim_defaults as sim_defaults
from halotools.sim_manager import CachedHaloCatalog
from halotools.empirical_models import PrebuiltHodModelFactory

# setup globals
halocat = None
HODmodel = None
catLbox = None
catboxz = None
# newLbox = None
cosmo = None
H0 = None
Om0 = None


def get_galtbl(getas='HOD'):
    """ getas='DF' returns the galaxy_table as a DF using galtbl_to_DF().
                            Only columns are galaxy phase space!
        else returns HODmodel.mock.galaxy_table (astropy table).
    """
    global HODmodel
    global catLbox
    global catboxz
    load_popmock() # (re)populate mock, set catLbox, catboxz
    galaxy_table = HODmodel.mock.galaxy_table
    if getas == 'DF':
        galaxy_table = galtbl_to_DF(galaxy_table)
    return [galaxy_table, catLbox, catboxz]



def galtbl_to_DF(galaxy_table, columns=None):
    """Returns galaxy_table (astropy Table) as DataFrame, keeping only columns."""
    if columns is None:
        columns = ['x','y','z', 'vx','vy','vz']
    df = galaxy_table.to_pandas()
    df_strip = df[columns]
    return df_strip



def load_popmock():
    """
    Loads a repopulated mock halo (need to set up default) to the system.
    Will load halocat and HODmodel if needed.
    """
    global halocat
    global HODmodel
    global catLbox
    global catboxz
    if HODmodel is None:
        print('\nSetting up halocat and HODmodel\n')
        halocat, HODmodel = setup_halosHOD() # fetch halo catalog and HODmodel, returns populated HODmodel.mock
        # HODmodel.mock.populate() # repopulate
        catLbox = halocat.Lbox[0]
        catboxz = halocat.redshift
    HODmodel.mock.populate() # repopulate



def load_cosmo(H0in=70.0, Om0in=0.3):
    """Load a global cosmology for consistency.
    """
    global cosmo
    global H0
    global Om0

    newcosmo = cosmology.FlatLambdaCDM(H0=H0in, Om0=Om0in)
    try: # pass if already loaded
        assert cosmo == newcosmo
    except: # set for first time, or with new params
        print('*** Loading global cosmo. \n\t Old cosmo = {}'.format(cosmo))
        H0 = H0in
        Om0 = Om0in
        cosmo = newcosmo
        print('\t New cosmo = {}'.format(cosmo))



# cosmo.comoving_distance(zred) returns "Comoving distance in Mpc to each input redshift.""
# from help(cosmo): Dimensionless Hubble constant: h = H_0 / 100 [km/sec/Mpc]
def push_box2z_old(galaxy_coords, redshift, Lbox):
    """Moves the x coordinates of galaxy_coords
        so the FACE of the box is at comoving_distance(redshift).
    Expects galaxy_coords as (gtot x 3) array with columns ['x','y','z'].
        extra columns at the end (e.g. velocities) are ok.
        all columns after the first one will be returned unchanged.
    Assumes coords are in Mpc/h
    Returns galaxy_coords with x coordinates shifted.
    """
    global cosmo
    if cosmo is None:
        load_cosmo()

    xx = galaxy_coords[:,0]
    # shift xx so that coordinates are strictly positive (i.e. move face to x=0)
    xx = xx+ Lbox/2
    # shift xx so the face is at comoving_distance(redshift)
    xzbox = (cosmo.comoving_distance(redshift).value)*cosmo.h # Mpc/h
    xx = xx+ xzbox

    # change given x coords and return the array (otherwise unchanged)
    galaxy_coords[:,0] = xx
    return galaxy_coords





def setup_halosHOD(simname='multidark', redshift=0.5, HODparams=None):
    """Loads a halo catalog from simname, redshift, and ROCKSTAR (hardcoded)
        (assumes catalog has been previously cached).

        Sets up a zheng07 (hardcoded) HOD model using HODparams
        (default values from Rongpu's LRG analysis (see Slack message)).

        Returns [halocat, HODmodel]
    """
    # Defaults
    if HODparams == None:
        HODparams=[12.89376701, 0.23939566, 12.26827089, 14.03372441, 1.32828278]
    sim_defaults.default_simname = simname
    sim_defaults.default_redshift = redshift # multidark actual z = 0.466
    halo_finder = 'rockstar'


    # Download halo catalogs (only needs to be done once)
    # from halotools.sim_manager import DownloadManager
    # dman = DownloadManager()
    # dman.download_processed_halo_table('multidark', 'rockstar', 0.5)



    # get Multidark + ROCKSTAR halo catalog
    halocat = CachedHaloCatalog(simname=simname, halo_finder=halo_finder, \
                version_name = 'halotools_v0p4', redshift=redshift)


    # zheng07 HOD model
    # https://halotools.readthedocs.io/en/latest/quickstart_and_tutorials/tutorials/model_building/preloaded_models/zheng07_composite_model.html
    HODmodel = PrebuiltHodModelFactory('zheng07', redshift=redshift) #, threshold = -21)
    HODmodel.param_dict['logMmin'] = HODparams[0] # Minimum mass required for a halo to host a central galaxy.
    HODmodel.param_dict['sigma_logM'] = HODparams[1] # Rate of transition from ⟨Ncen⟩=0⇒⟨Ncen=1⟩
    HODmodel.param_dict['logM0'] = HODparams[2] # Low-mass cutoff in ⟨Nsat⟩
    HODmodel.param_dict['logM1'] = HODparams[3] # Characteristic halo mass where ⟨Nsat⟩ begins to assume a power law form.
    HODmodel.param_dict['alpha'] = HODparams[4] # Power law slope of the relation between halo mass and ⟨Nsat⟩.
    HODmodel.param_dict



    # Create mock galaxy catalogs
    # (only last instance is available for future calls on HODmodel.mock)
    HODmodel.populate_mock(halocat) # use this command to create and populate the mock
    # HODmodel.mock.populate() # use this command to REpopulate

    return [halocat, HODmodel]
