""" Returns the inupt mock galaxy catalog with transformations:
        * Center box on origin then push the face to param_dict['z4push'].
        * Convert phase space coordinates to ra, dec, redshift.

    Use Args:
                param_dict = main.load_param_dict()
                gals_PS = gm.get_mock(param_dict, randoms=False)
"""

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from astropy.constants import c  # the speed of light



def transform(param_dict, gals_PS):
    """

    Returns
        gals_PS with transformed coordinates.
        gals_RDZ
        rands_PS
        rands_RDZ
        param_dict with zbin_edges and zbin_cens added.
    """

    p = param_dict.copy()

    gals_PS = stack_boxes_around_origin(p, gals_PS)
    gals_PS = shift_face_to_z4push(p, gals_PS)
    gals_RDZ, p['zbin_edges'], p['zbin_cens'] = get_ra_dec_z(p, gals_PS)

    return gals_PS, gals_RDZ, p






def stack_boxes_around_origin(param_dict, gals_PS):

    p = param_dict

    Nstack = p['mock_Nstack']
    assert ((Nstack%2 == 0) or (Nstack==1)), "param_dict['mock_Nstack'] should be 1 or even int."
    N2 = int(Nstack/2)

    if Nstack == 1: # shift origin to center and return
        gals_PS.loc[:,['x','y','z']] = gals_PS.loc[:,['x','y','z']] - p['mock_Lbox']/2.
        return gals_PS

    # get df of columns that will not be transformed
    static_cols_df = gals_PS[[c for c in list(gals_PS.columns.values) if c not in ['x','y','z']]]

    idx_omag = np.floor(np.log10(np.max(gals_PS.index.values))) # order of mag of max index
    nblist = [] # list to hold DataFrames for each new box
    nLin = np.linspace(-N2,N2-1,Nstack)*p['sim_Lbox']  # array of ints -N2 to N2-1
    deltax,deltay,deltaz = np.meshgrid(nLin,nLin,nLin) # holds deltas for each new box
    for b in range(Nstack**3): # for each new box
        boxdf = static_cols_df.copy(deep=True)

        # add index offset to ensure unique indices for each galaxy
        boxdf.index = boxdf.index + np.int32(b*10**(idx_omag+1))

        # create coords in new box
        dx,dy,dz = deltax.flat[b], deltay.flat[b], deltaz.flat[b]
        boxdf['x'] = gals_PS.x + dx
        boxdf['y'] = gals_PS.y + dy
        boxdf['z'] = gals_PS.z + dz

        nblist.append(boxdf)

    gals_PS = pd.concat(nblist, ignore_index=False)

    return gals_PS


def shift_face_to_z4push(param_dict, PSdf):
    """
    Args:
        PSdf  (df): Phase space coords dataframe.
                    Min cols 'x','y','z'.
                    Centered around origin.

    Returns:
        PSdf  (df): Input with 'x' column shifted to x(z4push)
    """

    p = param_dict
    z4push = p['sim_redshift'] if p['z4push']=='sim' else p['z4push']

    z_deltax = (p['cosmo'].comoving_distance(z4push).value) * p['cosmo'].h # Mpc/h
    PSdf.loc[:,'x'] = PSdf.loc[:,'x'] + p['mock_Lbox']/2. + z_deltax

    return PSdf


def get_ra_dec_z(param_dict, PSdf):
    """ PSdf = DataFrame with minimum columns {x,y,z,[h^-1 Mpc] vx,vy,vz}, rows=galaxies

    Returns:

        rdz           (DF): {RA, DEC, Redshift, zbin} with ra, dec in degrees
        zbin_edges (array): zbin bin edges
        zbin_cens  (array): zbin bin centers

    """
    p = param_dict
    c_km_s = c.to('km/s').value # speed of light

    # remove h scaling from position so we can use cosmo.comoving_distance
    x,y,z = ( PSdf.loc[:,ax] / p['cosmo'].h for ax in ['x','y','z'] )

    r = np.sqrt(x**2+y**2+z**2) # comoving distance from observer
    ct = z/r
    st = np.sqrt(1.0 - ct**2)
    cp = x/np.sqrt(x**2 + y**2)
    sp = y/np.sqrt(x**2 + y**2)

    ### RA, DEC
    # calculate spherical coordinates
    theta = np.arccos(z/r)
    phi = np.arctan2(y, x)
    # convert spherical coordinates into ra,dec
    dec = np.degrees(np.pi/2.0 - theta)
    ra = np.degrees(phi)
    # convert ra to interval [0,360] for calc_wtheta
    ra[ra<0.] = 360.+ ra

    ### Observed Redshift
    vr = PSdf.vx*st*cp + PSdf.vy*st*sp + PSdf.vz*ct # = radial peculiar comoving velocity
    # compute cosmological redshift and add contribution from peculiar velocity
    yy = np.arange(0, 3, 0.001)
    xx = p['cosmo'].comoving_distance(yy).value
    f = interp1d(xx, yy, kind='cubic')
    z_cos = f(r)
    redshift = z_cos+(vr/c_km_s)*(1.0+z_cos)
    # bin redshifts
    zbin_edges, zbin_cens = set_zbins(p, redshift) if 'zbin_edges' not in p.keys() \
                            else (p['zbin_edges'], p['zbin_cens'])
    zbin_loc = np.digitize(redshift, zbin_edges, right=True)
    zbin = np.array([ zbin_cens[i-1] for i in zbin_loc ])

    # Create DF with dtype float32.
    # This ensures that Corrfunc gets galaxies and randoms with same type (as required).
    # Use float32 since this is the output from HODmodel.mock.galaxy_table.
    rdz = pd.DataFrame(data={'RA':ra, 'DEC':dec, 'Redshift':redshift, 'zbin':zbin},
                        dtype='float32')

    return rdz, zbin_edges, zbin_cens


def set_zbins(param_dict, redshifts):
    """ Generates redshift bins based on input redshifts.
        Returns bin edges and bin centers.
    """
    zmin, zmax = np.floor(10*redshifts.min())/10, np.ceil(10*redshifts.max())/10
    w = param_dict['zbin_width']
    num_zbins = np.ceil((zmax - zmin)/ w)

    edges = np.array([ zmin+ w*i for i in range(int(num_zbins)+1) ])
    cens = np.array([(edges[i]+edges[i+1])/2 for i in range(len(edges)-1)])

    return edges, cens
