import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

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




# NOT WORKING YET
# look at galaxy distribution
def plot_galaxies(galaxy_table, gal_frac=0.05):
    # get a random sample
    lg = np.arange(len(galaxy_table['x']))
    np.random.shuffle(lg)
    lg = lg[:int(0.05*len(galaxy_table['x']))]

    plt.figure(figsize=(8,8))
    ax = plt.axes(projection='3d')
    ax.scatter3D(np.asarray(galaxy_table['x'][lg]), np.asarray(galaxy_table['y'][lg]), \
                    np.asarray(galaxy_table['z'][lg]), s=0.1)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()
