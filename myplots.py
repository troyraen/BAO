import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import imp

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
    plt.show()




# look at galaxy distribution
def plot_galaxies(galaxy_table, gal_frac=0.05, coords='xyz'):
    """ Plots 3D galaxy distribution using coords cordinate basis.
            if coords == 'radecz', will call get_ra_dec_z to transform the coordinates
    """
    # get a random sample
    lg = np.arange(len(galaxy_table['x']))
    np.random.shuffle(lg)
    lg = lg[:int(0.05*len(galaxy_table['x']))]

    plt.figure(figsize=(10,6))
    ax = plt.axes(projection='3d')

    if coords == 'xyz':
        x = np.asarray(galaxy_table['x'][lg])
        y = np.asarray(galaxy_table['y'][lg])
        z = np.asarray(galaxy_table['z'][lg])
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

    elif coords == 'radecz':
        x, y, z = get_ra_dec_z(galaxy_table)
        ax.set_xlabel('RA')
        ax.set_ylabel('DEC')
        ax.set_zlabel('z')

    else:
        raise Exception('coords must be either \'xyz\' or \'radecz\'\n\t {} not a recognized option'.format(coords))

    ax.scatter3D(x, y, z, s=0.1)
    plt.show()
