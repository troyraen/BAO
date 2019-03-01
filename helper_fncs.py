import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from astropy import cosmology
from astropy.constants import c  # the speed of light

import setup as su




# SDSS DR10 max photo zErr = 0.365106
def bin_redshifs(rdz, zspace = 0.365, validate=False):
    """rdz = DataFrame with minimum columns {RA, DEC, Redshift} (as returned by get_ra_dec_z)
        zspace = desired spacing between bins (result will not be exact)

    Returns:
        Original DataFrame with added column
            containing the bin (zbcens value) the galaxy resides in.
        Array of bin edges
    """
    zmin = rdz.Redshift.min()
    zmax = rdz.Redshift.max()
    # eps = 0.01
    # zbins = np.arange(zmin, zmax+eps, zspace)
    num_binedges = int(np.ceil((zmax-zmin)/zspace))
    zbins = np.linspace(zmin,zmax,num_binedges)
    zbcens = (zbins[:-1]+zbins[1:])/2
    # create zbin masks for rdz dataframe
        # add a column to rdz containing the zbin (zbcens value) the galaxy resides in
    # given z, find which bin it's in and get the value of the bin center
    rdz['zbin'] = rdz['Redshift'].apply(find_bin_center, **{"bin_edges":zbins, "bin_centers":zbcens})
    # rdz.apply(lambda inval: hf.find_bin_center(inval, zbins, zbcens), rdz.Redshift)

    # make sure the operation worked as expected:
    if validate:
        for index, row in rdz.iterrows():
            rdzbin = row.zbin
            truezbin = find_bin_center(row.Redshift, bin_edges=zbins, bin_centers=zbcens)
            assert(rdzbin==truezbin)


    return [rdz, zbins]



def find_bin_center(inval, bin_edges=None, bin_centers=None):
    """Returns the bin_centers value corresponding to the
        bin (defined by bin_edges) that inval (scalar) falls in to.
    """
    for i in range(len(bin_centers)):
        if (bin_edges[i] <= inval) and (inval < bin_edges[i+1]):
            return bin_centers[i]
    if inval == bin_edges[-1]: # inval == top edge of last bin
        return bin_centers[-1]
    # if you get here, inval is not in any of the bins
    print('{} did not fall within any bins.'.format(inval))
    return None

def get_ra_dec_z(galdf, usevel=True):
    """Most of this is taken from Duncan Campbell's function mock_survey.ra_dec_z
        galdf should be a DataFrame with minimum columns {x,y,z, vx,vy,vz}
        usevel = True will add reshift due to perculiar velocities
        Returns galdf columns {RA, DEC, Redshift} added (or updated) with ra, dec in degrees
    """

    # Get new DF with RA,DEC,Z. galdf indexing is preserved.
    rdz = galdf.apply(get_ra_dec_z_calculate, axis=1)
    # Join it to the original galdf and return it
    galdf = galdf.join(rdz)
    return galdf

# First, interp redshift once
yy = np.arange(0, 2.0, 0.001)
su.load_cosmo()
xx = su.cosmo.comoving_distance(yy).value
f = interp1d(xx, yy, kind='cubic')
def get_ra_dec_z_calculate(gal, usevel=True):
    """ This is intended to be called using galdf.apply(get_ra_dec_z_calculate, axis=1)
        Returns a Series with indices {RA,DEC,Redshift} for the given row.
    """


    x = np.array([gal.x, gal.y, gal.z])
    v = np.array([gal.vx, gal.vy, gal.vz])

# calculate the observed redshift
    c_km_s = c.to('km/s').value

    # remove h scaling from position so we can use the cosmo object
    x = x/su.cosmo.h

    # compute comoving distance from observer
    r = np.sqrt(x[0]**2+x[1]**2+x[2]**2)

    # compute radial velocity
    ct = x[2]/r
    st = np.sqrt(1.0 - ct**2)
    cp = x[0]/np.sqrt(x[0]**2 + x[1]**2)
    sp = x[1]/np.sqrt(x[0]**2 + x[1]**2)
    vr = v[0]*st*cp + v[1]*st*sp + v[2]*ct # = radial peculiar velocity (comoving)?

    # compute cosmological redshift and add contribution from perculiar velocity
    # yy = np.arange(0, 2.0, 0.001)
    # xx = su.cosmo.comoving_distance(yy).value
    # f = interp1d(xx, yy, kind='cubic')
    z_cos = f(r)
    redshift = z_cos+(vr/c_km_s)*(1.0+z_cos) if usevel else z_cos

    # calculate spherical coordinates
    theta = np.arccos(x[2]/r)
    phi = np.arctan2(x[1], x[0])

    # convert spherical coordinates into ra,dec
    dec = np.degrees(np.pi/2.0 - theta)
    ra = np.degrees(phi)
    # convert ra to interval [0,360] for calc_wtheta
    ra = 360+ ra if ra<0 else ra
    # print(max(ra), min(ra))
    # msk = (np.ma.masked_less(ra,0)).mask # True for all values of ra < 0
    # ra[msk] = 360+ ra[msk]
    # print(max(ra), min(ra))

    # collect results
    # rdz = np.vstack((ra,dec,redshift)).T
    # rdz = pd.DataFrame(rdz, columns=['RA','DEC','Redshift'])

    # coords = np.vstack([galaxy_table['x'], galaxy_table['y'], galaxy_table['z']]).T # check these units, mock_survey.ra_dec_z expects Mpc/h
    # vels = np.vstack([galaxy_table['vx'], galaxy_table['vy'], galaxy_table['vz']]).T # mock_survey.ra_dec_z expects km/s
    #
    # ra, dec, z = mock_survey.ra_dec_z(coords, vels) # returns ra, dec in radians
    # ra = np.degrees(ra) # [0, 90] degrees
    # dec = np.degrees(dec) # [-90, 0] degrees
    #
    # return [ra, dec, z]

    rdzdic = {'RA':ra, 'DEC':dec, 'Redshift':redshift}
    rdz = pd.Series(data=rdzdic)
    return rdz




def get_ra_dec_z_old(ps_coords, usevel=True):
    """Most of this is taken from Duncan Campbell's function mock_survey.ra_dec_z
        ps_coords should be ndarray (ngals x 6) (columns = {x,y,z, vx,vy,vz})
        usevel = True will add reshift due to perculiar velocities
        Returns dataframe ngals x 3 {RA, DEC, Redshift} with ra, dec in degrees
    """

    x = ps_coords[:,0:3]
    v = ps_coords[:,3:6]

# calculate the observed redshift
    c_km_s = c.to('km/s').value

    # remove h scaling from position so we can use the cosmo object
    x = x/su.cosmo.h

    # compute comoving distance from observer
    r = np.sqrt(x[:, 0]**2+x[:, 1]**2+x[:, 2]**2)

    # compute radial velocity
    ct = x[:, 2]/r
    st = np.sqrt(1.0 - ct**2)
    cp = x[:, 0]/np.sqrt(x[:, 0]**2 + x[:, 1]**2)
    sp = x[:, 1]/np.sqrt(x[:, 0]**2 + x[:, 1]**2)
    vr = v[:, 0]*st*cp + v[:, 1]*st*sp + v[:, 2]*ct # = radial peculiar velocity (comoving)?

    # compute cosmological redshift and add contribution from perculiar velocity
    yy = np.arange(0, 4.0, 0.001)
    xx = su.cosmo.comoving_distance(yy).value
    f = interp1d(xx, yy, kind='cubic')
    z_cos = f(r)
    redshift = z_cos+(vr/c_km_s)*(1.0+z_cos) if usevel else z_cos

    # calculate spherical coordinates
    theta = np.arccos(x[:, 2]/r)
    phi = np.arctan2(x[:, 1], x[:, 0])

    # convert spherical coordinates into ra,dec
    dec = np.degrees(np.pi/2.0 - theta)
    ra = np.degrees(phi)
    # convert ra to interval [0,360] for calc_wtheta
    # print(max(ra), min(ra))
    msk = (np.ma.masked_less(ra,0)).mask # True for all values of ra < 0
    ra[msk] = 360+ ra[msk]
    # print(max(ra), min(ra))

    # collect results
    rdz = np.vstack((ra,dec,redshift)).T
    rdz = pd.DataFrame(rdz, columns=['RA','DEC','Redshift'])

    # coords = np.vstack([galaxy_table['x'], galaxy_table['y'], galaxy_table['z']]).T # check these units, mock_survey.ra_dec_z expects Mpc/h
    # vels = np.vstack([galaxy_table['vx'], galaxy_table['vy'], galaxy_table['vz']]).T # mock_survey.ra_dec_z expects km/s
    #
    # ra, dec, z = mock_survey.ra_dec_z(coords, vels) # returns ra, dec in radians
    # ra = np.degrees(ra) # [0, 90] degrees
    # dec = np.degrees(dec) # [-90, 0] degrees
    #
    # return [ra, dec, z]

    return rdz
