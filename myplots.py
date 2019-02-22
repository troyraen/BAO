import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import imp

from astropy import cosmology
import calc_wtheta as cw

# set plot defaults
mpl.rcParams['font.size'] = 14
mpl.rcParams['legend.fontsize'] = 'small'
mpl.rcParams['figure.titlesize'] = 'medium'



def plot_wtheta(bcens, wtheta):
    plt.figure()
    plt.scatter(bcens, wtheta)
    plt.plot(bcens, wtheta)
    # plt.xlim(2.5,20)
    # plt.ylim(-0.0015, 0.005)
    # plt.legend()
    plt.xlabel(r'$\theta$ [deg]')
    plt.ylabel(r'w($\theta$)')
    plt.title(r'w($\theta$)')
    plt.tight_layout()
    # plt.savefig('./wtheta.png')
    plt.show(block=False)




# look at galaxy distribution
def plot_galaxies(galaxy_table, gal_frac=0.05, coords='xyz'):
    """ Plots 3D galaxy distribution using coords cordinate basis.
            if coords == 'radecz', will call get_ra_dec_z to transform the coordinates
    """
    # get a random sample
    lg = np.arange(len(galaxy_table['x']))
    np.random.shuffle(lg)
    lg = lg[:int(gal_frac*len(galaxy_table['x']))]

    plt.figure(figsize=(13,13))
    ax = plt.axes(projection='3d')

    if coords == 'xyz':
        x = np.asarray(galaxy_table['x'][lg])
        y = np.asarray(galaxy_table['y'][lg])
        z = np.asarray(galaxy_table['z'][lg])
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

    elif coords == 'radecz':
        x, y, z = cw.get_ra_dec_z(galaxy_table)
        x = x[lg]
        y = y[lg]
        z = z[lg]
        ax.set_xlabel('RA')
        ax.set_ylabel('DEC')
        ax.set_zlabel('z')
        ax.view_init(azim=260, elev=95) # rotate the view to physical line of sight

    else:
        raise Exception('coords must be either \'xyz\' or \'radecz\'\n\t {} not a recognized option'.format(coords))

    ax.scatter3D(x, y, z, s=1)
    plt.show(block=False)



# plot comoving distance vs redshift
# call with:
    # halocat, HODmodel = sm.setup_halosHOD() # fetch halo catalog and HODmodel, returns populated HODmodel.mock
    # zred = halocat.redshift
    # plot_dist_redshift(log=False, fout=False, zred=zred)
def plot_dist_redshift(zspace=None, log=False, fout=None, zred=None):
    """ pass zspace array for plotting
        pass zred = redshift of box to put point on plot
        fout = path as string to save file
    """
    cosmo = cosmology.FlatLambdaCDM(H0=70.0, Om0=0.3)
    zz = zspace if (zspace is not None) else np.arange(0.01, 2.0, 0.001)
    xx = cosmo.comoving_distance(zz).value # Mpc
    xxh = cosmo.comoving_distance(zz).value*cosmo.h # Mpc/h
    plt.figure()
    if log:
        plt.loglog(zz, xx, label='Mpc (h=0.7)')
        plt.loglog(zz, xxh, label='Mpc/h')
    else:
        plt.plot(zz, xx, label='Mpc (h=0.7)')
        plt.plot(zz, xxh, label='Mpc/h')
    if zred is not None:
        xzbox = (cosmo.comoving_distance(zred).value)*cosmo.h # Mpc/h
        plt.scatter(zred, xzbox, label='zbox= {}'.format(zred))
    plt.xlabel('Redshift')
    plt.ylabel('Comoving Distance')
    plt.legend()
    plt.grid(linestyle='-', linewidth='0.5', color='0.7')
    if fout is not None:
        plt.savefig(fout)
    plt.show(block=False)
