import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
# import imp

from astropy import cosmology
import calc_wtheta as cw

# set plot defaults
mpl.rcParams['figure.figsize'] = [8.0, 6.0]
mpl.rcParams['font.size'] = 14
mpl.rcParams['legend.fontsize'] = 'small'
mpl.rcParams['figure.titlesize'] = 'medium'



def getplot_zruntimes(zrunfout='data/zruntime.dat'):
    # get and plot zrun calc times
    zrf = pd.read_csv(zrunfout, delim_whitespace=True)
    # size = np.log10(zrf.index.values+1.1)*1000
    size = zrf.index.values*50
    zrf.plot.scatter(x='nthreads',y='calctime', s=size, \
            c='numgals', colormap='YlGnBu', alpha=.5 )
    # plt.semilogy()
    plt.title('zbin calc_wtheta runtimes.')
    plt.tight_layout()
    plt.show(block=False)
    return zrf




def plot_wtheta(wdf):
    """wdf = DataFrame with columns wtheta(theta_bcens), 'zbin', 'mock'
        (if multiple mock nums for given zbin, get average.)
        Assumes all column names that can be converted to a float are theta bin centers
        Plots wtheta(theta_bcens), one line for each zbin
    """

    bincols = cw.get_tbins(wdf) # get list of theta bin column names
    # create new df of wtheta values using pivot_table
        # calculating mean wtheta for each zbin.
        # Transpose to use df.plot()
    wdf.zbin = np.round(wdf.zbin,1) # zbin centers are not exactly the same. this needs to be fixed
    wtheta = (pd.pivot_table(wdf[bincols+['zbin']], index='zbin')).T
    wtheta.rename(index=lambda c: np.double(c), inplace=True) # change index dtype to double

    plt.figure()
    wtheta.plot()
    plt.xlabel(r'$\theta$ [deg]')
    plt.ylabel(r'$w(\theta)$')
    plt.tight_layout()
    plt.show(block=False)
    # for idx, row in wtheta.iterrows():
    #     print(row[bincols[0]])
    #     # ax.scatter(row[1], loc, label=row[1].name)

    # plt.figure()
    # tbins = np.asarray(bincols,dtype=np.double)
    # zbcens = wdf.zbin.unique()
    # for zzz in zbcens:
    #     # if there are multiple mock nums for given zbin, get average
    #     wtheta = wdf[bincols].loc[wdf.zbin == zzz] # get wtheta columns of row corresponding to zzz
    #     plt.plot(tbins, wtheta, label='zbin = {0:1.5f}'.format(zzz))
    #
    # plt.xlabel(r'$\theta$ [deg]')
    # plt.ylabel(r'w($\theta$)')
    # plt.legend()
    # plt.title(r'w($\theta$)')
    # plt.tight_layout()
    # plt.show(block=False)




def plot_wtheta_old(bcens, wtheta):
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
def plot_galaxies(galaxies, gal_frac=0.05, title='Galaxies', coords='xyz'):
    """ Plots 3D galaxy distribution.
        galaxies assumed to be DataFrame
    """
    # get a random sample
    gs = galaxies.sample(int(len(galaxies)*0.05))
    # lg = np.arange(len(galaxies['x']))
    # np.random.shuffle(lg)
    # lg = lg[:int(gal_frac*len(galaxies['x']))]

    plt.figure()
    ax = plt.axes(projection='3d')
    # plt.figure().gca(projection='3d')

    ax.scatter3D(gs.x, gs.y, gs.z, s=1)
    plt.title(title)
    plt.tight_layout()
    plt.show(block=False)



# look at galaxy distribution
def plot_galaxies_old(galaxy_table, gal_frac=0.05, coords='xyz', title='Galaxies'):
    """ Plots 3D galaxy distribution using coords cordinate basis.
        galaxy_table assumed to be astropy table with cols {'x','y','z'}
        coords should be one of 'xyz' (h^-1 Mpc assumed),
                                NO!: 'xyred' (h^-1 Mpc assumed for x and y),
                                NO!:    'rdz' (degrees assumed for ra and dec)
    """
    # get a random sample
    lg = np.arange(len(galaxy_table['x']))
    np.random.shuffle(lg)
    lg = lg[:int(gal_frac*len(galaxy_table['x']))]

    plt.figure(figsize=(8,8))
    ax = plt.axes(projection='3d')

    if coords == 'xyz':
        x = np.asarray(galaxy_table['x'][lg])
        y = np.asarray(galaxy_table['y'][lg])
        z = np.asarray(galaxy_table['z'][lg])
        ax.set_xlabel('x h^-1 [Mpc]')
        ax.set_ylabel('y h^-1 [Mpc]')
        ax.set_zlabel('z h^-1 [Mpc]')
    # elif coords == 'xyred':
    #     ax.set_xlabel('x h^-1 [Mpc]')
    #     ax.set_ylabel('y h^-1 [Mpc]')
    #     ax.set_zlabel('redshift')
    # elif coords == 'rdz':
    #     ax.set_xlabel('x h^-1 [Mpc]')
    #     ax.set_ylabel('y h^-1 [Mpc]')
    #     ax.set_zlabel('redshift')


    # elif coords == 'radecz':
    #     x, y, z = cw.get_ra_dec_z(galaxy_table)
    #     x = x[lg]
    #     y = y[lg]
    #     z = z[lg]
    #     ax.set_xlabel('RA')
    #     ax.set_ylabel('DEC')
    #     ax.set_zlabel('z')
    #     ax.view_init(azim=260, elev=95) # rotate the view to physical line of sight

    else:
        raise Exception('coords must be either \'xyz\' or \'xyred\'\n\t {} not a recognized option'.format(coords))

    ax.scatter3D(x, y, z, s=1)
    plt.title(title)
    plt.tight_layout()
    plt.show(block=False)




# plot comoving distance vs redshift
# call with:
    # halocat, HODmodel = sm.setup_halosHOD() # fetch halo catalog and HODmodel, returns populated HODmodel.mock
    # zred = halocat.redshift
    # plot_dist_redshift(log=False, fout=False, zred=zred)
def plot_codist_redshift(zspace=None, log=False, fout=None, zred=None, cosmo=None):
    """ pass zspace array for plotting
        pass zred = redshift of box to put point on plot
        fout = path as string to save file
    """
    if cosmo is None:
        cosmo = cosmology.FlatLambdaCDM(H0=70.0, Om0=0.3)
    zz = zspace if (zspace is not None) else np.arange(0.01, 2.0, 0.001)
    xx = cosmo.comoving_distance(zz).value # Mpc
    xxh = cosmo.comoving_distance(zz).value*cosmo.h # Mpc/h
    plt.figure()
    if log:
        plt.loglog(zz, xx, label='Mpc (h=0.7)')
        plt.loglog(zz, xxh, label='Mpc/h')
    else:
        plt.plot(zz, xx, label='Mpc (h={})'.format(cosmo.h))
        plt.plot(zz, xxh, label='Mpc/h')
    if zred is not None:
        xzbox = (cosmo.comoving_distance(zred).value)*cosmo.h # Mpc/h
        plt.scatter(zred, xzbox, label='zbox= {}'.format(zred))
    plt.xlabel('Redshift')
    plt.ylabel('Comoving Distance')
    plt.legend()
    plt.grid(linestyle='-', linewidth='0.5', color='0.7')
    plt.tight_layout()
    if fout is not None:
        plt.savefig(fout)
    plt.show(block=False)
