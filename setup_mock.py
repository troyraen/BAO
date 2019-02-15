import numpy as np

import halotools.sim_manager.sim_defaults as sim_defaults
from halotools.sim_manager import CachedHaloCatalog
from halotools.empirical_models import PrebuiltHodModelFactory


def stack_boxes(galaxy_table, Nstack=20, Lbox=1000):
    """Takes object HODmodel.mock.galaxy_table (gtin, assumed coordinates {x,y,z} in [0,Lbox]) and
            stacks Nstack (must be an even integer) of them together
        Returns a np.array (gtot x 6) with {x,y,z} generated periodically, {vx,vy,vz} same as original point
    """

    N2 = int(Nstack/2)
    assert Nstack%2 == 0, 'Nstack should be an even integer.'

    gt = galaxy_table
    x = gt['x']; y = gt['y']; z = gt['z']
    vx = gt['vx']; vy = gt['vy']; vz = gt['vz']

    ng = len(x) # number of galaxies in original box
    gtot = ng* Nstack**3 # new total # galaxies
    newgals = np.zeros((gtot,6))
    nLin = np.linspace(-N2,N2-1,Nstack)*Lbox  # array of ints -N2 to N2-1
    # boxes will be stacked around the origin by adding nLin to each galaxy coordinate
    gnew = 0
    for g in range(ng): # for each galaxy
        vg = [ vx[g], vy[g], vz[g] ] # get x,y,z velocities (same for each created coord)
        # create coordinates for stacked boxes in each direction
        xstk = x[g]*np.ones(Nstack) +nLin # array of length Nstack
        ystk = y[g]*np.ones(Nstack) +nLin
        zstk = z[g]*np.ones(Nstack) +nLin
        for i in xstk:
            for j in ystk:
                for k in zstk:
                    coords = np.array([ i, j, k ] + vg)
                    newgals[gnew] = coords # add the new galaxy
                    gnew = gnew+1

    return newgals




def get_ra_dec_z(galaxy_table):
    coords = np.vstack([galaxy_table['x'], galaxy_table['y'], galaxy_table['z']]).T # check these units, mock_survey.ra_dec_z expects Mpc/h
    vels = np.vstack([galaxy_table['vx'], galaxy_table['vy'], galaxy_table['vz']]).T # mock_survey.ra_dec_z expects km/s

    ra, dec, z = mock_survey.ra_dec_z(coords, vels) # returns ra, dec in radians
    ra = np.degrees(ra) # [0, 90] degrees
    dec = np.degrees(dec) # [-90, 0] degrees

    return [ra, dec, z]



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
