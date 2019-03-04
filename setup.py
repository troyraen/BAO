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
newLbox = None
cosmo = None
H0 = None
Om0 = None


def get_galtbl(getas='HOD'):
    """ getas='DF' returns the galaxy_table as a DF using galtbl_to_DF().
                            Only columns are galaxy phase space!
        else returns HODmodel.mock.galaxy_table (astropy table).
    """
    global HODmodel
    if HODmodel is None:
        load_popmock()
    galaxy_table = HODmodel.mock.galaxy_table
    if getas == 'DF':
        galaxy_table = galtbl_to_DF(galaxy_table)
    return galaxy_table



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

    # try:
    #     halocat and HODmodel
    # except:
    #     print('\nSetting up halocat and HODmodel\n')
    #     halocat, HODmodel = setup_halosHOD() # fetch halo catalog and HODmodel, returns populated HODmodel.mock
    #     # HODmodel.mock.populate() # repopulate
    #     catLbox = halocat.Lbox[0]
    #     catboxz = halocat.redshift

    # for nthreads in [48, 32, 12]:
    # print('\nDoing nthreads = {}\n'.format(nthreads))
    HODmodel.mock.populate() # repopulate



def load_cosmo(H0in=70.0, Om0in=0.3):
    """Load a global cosmology for consistency.
    """
    global cosmo
    global H0
    global Om0
    H0 = H0in
    Om0 = Om0in

    cosmo = cosmology.FlatLambdaCDM(H0=H0, Om0=Om0)
    print('*** cosmo changed to {}'.format(cosmo))


# cosmo.comoving_distance(zred) returns "Comoving distance in Mpc to each input redshift.""
# from help(cosmo): Dimensionless Hubble constant: h = H_0 / 100 [km/sec/Mpc]
def push_box2z(galdf, redshift, Lbox):
    """ Moves the box so the x-FACE is at comoving_distance(redshift).
        Returns galdf with only 'x' column changed.
    """
    global cosmo
    if cosmo is None:
        load_cosmo()

    # Shift x so coordinates are strictly positive (i.e. move face to x=0)
    # Then push the face to comoving_distance(redshift)
    deltax = Lbox/2. + (cosmo.comoving_distance(redshift).value)*cosmo.h # Mpc/h
    galdf.x = galdf.x + deltax

    return galdf


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



def stack_boxes_gen_new_rows(gal, Nstack=None, ogLbox=None):
    """ This gets applied to each row of the original dataframe.
        gal = single galaxy data as series (single row from original df)
        Generates new 'xyz' coordinates by stacking Nstack**3 boxes.
        Copy all other columns.
        Returns a DataFrame with Nstack**3 galaxies (rows)
    """
    N2 = int(Nstack/2)
    assert Nstack%2 == 0, 'Nstack should be an even integer.'

    dfcols = list(gal.index.values)
    extra_vals = gal[[c for c in dfcols if c not in ['x','y','z']]] # get series of just extra values

    numcols = len(dfcols)
    gtot = Nstack**3 # num galaxies generated from single original galaxy
    newgals = np.zeros((gtot,numcols)) # create an array to hold the output
    nLin = np.linspace(-N2,N2-1,Nstack)*ogLbox  # array of ints -N2 to N2-1
    # boxes will be stacked around the origin by adding nLin to each galaxy coordinate
    # create coordinates for stacked boxes in each direction
    xstk = gal.x*np.ones(Nstack) +nLin # array of length Nstack
    ystk = gal.y*np.ones(Nstack) +nLin
    zstk = gal.z*np.ones(Nstack) +nLin
    gnew = 0
    for i in xstk:
        for j in ystk:
            for k in zstk:
                coords = np.array([ i, j, k ] + list(extra_vals))
                newgals[gnew] = coords # add the new galaxy
                gnew = gnew+1
    ngdf = pd.DataFrame(newgals, columns=dfcols)
    return ngdf



def stack_boxes(galdf, Nstack=2, ogLbox=1000):
    """
    galdf = DataFrame representing a mock box,
            with minimum columns {'x','y','z', 'vx','vy','vz'},
            (assumes coordinates {x,y,z} in [0,Lbox]).
    Stacks Nstack^3 (must be an EVEN integer) galdf mock boxes together around the origin.
    Sets global newLbox.
    Returns a DataFrame of galaxies with the origin at the center of the box and boxsize=newLbox.
            columns {x,y,z} generated periodically, {vx,vy,vz, all others} same as original galaxies.
    """
    # Setup
    global newLbox
    # global catLbox
    newLbox = ogLbox*Nstack
    print('*** Warning: su.stack_boxes assumes the original box is strictly in the 1st quadrant with the origin at the corner. ***')

    if Nstack == 0: # Just move the origin to the center of the box.
        print('Moving origin to box center...')
        newLbox = ogLbox
        L2 = newLbox/2.
        newgals = galdf
        newgals.x = newgals.x-L2; newgals.y = newgals.y-L2; newgals.z = newgals.z-L2

    else:
        # Iterate over rows of galdf.
        # Generate new 'xyz' coordinates. Copy all other columns.
        nglist = [] # list to hold each new DataFrame generated from a single galaxy
        for idx, gal in galdf.iterrows():
            nglist.append(stack_boxes_gen_new_rows(gal, Nstack=Nstack, ogLbox=ogLbox))
        # Create new df from dfs in nglist
        newgals = pd.concat(nglist, ignore_index=True)

    return newgals




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
