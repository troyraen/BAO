import numpy as np
import pandas as pd
import time
import datetime
from scipy.interpolate import interp1d
from pathlib import Path
import os

from astropy import cosmology
from astropy.constants import c  # the speed of light

import setup as su


def file_ow(fin):
    """ Moves file fin.txt to fin_ow_%m%d%y_%H%M.txt
        Returns the new file name.
    """

    # split file name
    fsplit = fin.split('.')
    beg, end = '.'.join(fsplit[0:-1]), '.'+fsplit[-1]

    # get date and time to use in new file name
    dtm = datetime.datetime.now()
    owdtm = dtm.strftime("_ow_%m%d%y_%H%M")

    # create new file name
    fout = beg + owdtm + end

    # move the file
    print('Moving existing file {0} to {1}'.format(fin, fout))
    os.rename(fin, fout)

    return fout



def write_report_times(report_times, fname):
    """ report_times should be a dict with key = column_name, val = column_value
        Writes values to fname.
            Moves an existing file if columns don't match.
    """
    # Setup
    rt = report_times
    fcol_width = rt.pop('fcol_width', 20) # = 20 if fcol_width not in ft.keys
    rtcols = list(rt.keys())
    rtvals = list(rt.values())
    # print(rtvals)

    fpath = Path(fname)
    if fpath.is_file(): # Check the header
        fcols = list(np.genfromtxt(fname, max_rows=1, dtype=str))
        try:
            assert rtcols == fcols
        except: # columns don't match. Move the file
            print('*** report_times columns do not match existing file {}. ***'.format(fname))
            file_ow(fname)

    # If fname has been moved (above) or never existed, create new file.
    if not fpath.is_file():
        hdr = ''.join(str(x).rjust(fcol_width) for x in rtcols)
        print(hdr, file=open(fname, 'w'))

    # Now fname exists, so append data
    rtdat = ''.join(str(x).rjust(fcol_width) for x in rtvals)
    print(rtdat, file=open(fname, 'a'))

    return


def time_code(start, unit=None):

    """ Usage:  Use an OrderedDict, rt, to collect code runtimes as follows.
                Then use write_report_times(rt, file_name) to write everything to a file.

rt['CODE_NAME'] = hf.time_code('start') #.TS. get code start time
"<<-- Code you want to time goes here. -->> " #.TC.
rt['CODE_NAME'] = hf.time_code(rt['CODE_NAME'], unit='min') #.TE. replace start time with runtime in minutes

        start == 'start' returns system time as a float.
        Pass the result of time_code(start) as the argument in the next call to the function.
        type(start) == float => returns the runtime in units {'sec', 'min', 'hr', 'day'}.
                                            Assumes start is the start time.
    """
    if start == 'start':
        return time.time() # float

    elif type(start) == float:
        stop = time.time() # time the function
        runtime = (stop-start) # in seconds

        if unit == None or unit == 'min':
            return np.round(runtime/60., 2) # in minutes
        elif unit == 'sec':
            return runtime # in seconds
        elif unit == 'hr':
            return runtime/3600. # in hours
        elif unit == 'day':
            return runtime/3600./24. # in days
        else:
            print('\n*** time_code given invaid unit. (start={0}, unit={1}.) Returning minutes. ***\n'.format(start,unit))
            return runtime/60. # in minutes

    else:
        print('\n*** time_code received invalid argument, start = {}. ***\n\t*** Returning 0. ***\n'.format(start))
        return 0


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
        usevel = True will add reshift due to peculiar velocities
        Returns DF with columns {RA, DEC, Redshift} with ra, dec in degrees
    """

    # Get new DF with RA,DEC,Z. galdf indexing is preserved.
    rdz = galdf.apply(get_ra_dec_z_calculate, axis=1, **{'usevel':usevel})
    # Join it to the original galdf and return it
    # galdf = galdf.join(rdz)
    return rdz


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
