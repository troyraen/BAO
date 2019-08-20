import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def plot_stats(fdat, save=None, show=True):
    """ Plots stats (1 per column) in file fdat.
        Args:
        fdat (string): path to stats.dat file as written by MockBox.write_stat_to_file()
    """
    holdfsz = mpl.rcParams['figure.figsize'] # keep to reset later
    mpl.rcParams['figure.figsize'] = [14.0, 4.0]

    df = load_statsdat(fdat)

    sdf = df.groupby('statname').mean() # df
    sdf['NR/NG'] = (sdf['Nrands']/sdf['Ngals']).astype(int)
    # validate_statmeans_xbins(df) # make sure we haven't averaged different xbins
    lendf = df.groupby('statname').size() # series with # of mocks aggragated in each df above

    nrows, ncols = 1, len(lendf)
    fig, axs = plt.subplots(nrows, ncols, sharex=False, sharey=False)

    # Plot for each stat
    for i, (stat, row) in enumerate(sdf.iterrows()):
        ax = axs[i]

        x,y,ylabel = get_bins_stats(row, stat)
        try:
            lbl = r'{:.2f}$\pm${:.2f}, {:.0f}, {}, {}'.format(\
                        row.zbin, row.zwidth/2, row.Nstack, lendf.loc[stat], row['NR/NG'])
        except: # use for old format files with no zwidth
            lbl = r'{:.2f}, {:.0f}, {}, {}'.format(\
                        row.zbin, row.Nstack, lendf.loc[stat], row['NR/NG'])
        ax.semilogx(x,y, label=lbl)

        ax.axhline(0, c='0.5')
        ax.grid(which='both')
        ax.set_title(stat)
        ax.legend(title='zbin, Nstk, Nm, NR')
        ax.set_xlabel(r'$\theta$ [deg]' if stat not in ['wp','xi'] else r'r $h^{-1}$ [Mpc]')
        ax.set_ylabel(ylabel)

        # if stat != 'wtheta':
        #     ax.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
        #     ax.xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))

    plt.tight_layout()
    if save is not None:
        plt.savefig(save)
    if show:
        plt.show(block=False)

    mpl.rcParams['figure.figsize'] = holdfsz
    return None


def load_statsdat(fdat, stat=None, clean=False):
    """ Load correlation stats data from file fdat, as written by MockBox.write_stat_to_file.
        Returns DataFrame of file data.

        fdat (string or list of strings):   stats output file path(s)
        stat (string):  value in column statname in fdat
        clean (bool):   == True will drop rows with Nrands<=1000
    """

    if type(fdat)==str:
        df = pd.read_csv(fdat, delim_whitespace=True, comment='#')
    elif type(fdat)==list:
        dflst = []
        for f in fdat:
            dflst.append(pd.read_csv(f, delim_whitespace=True, comment='#'))
        df = pd.concat(dflst)
    else:
        print('mp.load_statsdat() received invalid argument.')
        print('fdat should be path(s) to stats files as string or list of strings')
        return None

    if stat is not None:
        df = df.query("statname=='{}'".format(stat))
    if clean:
        df = df.query("Nrands>1000")
    return df


def validate_statmeans_xbins(df, avg_zbins=False):
    """ Checks that xbins of collapsed statname have std_dev == 0
            and therefore are all the same.
        avg_zbins == True allows for each zbin to have unique theta bins
    """

    errmsg = "Operation groupby 'statname' has averaged values in different theta or r bins."
    if avg_zbins:
        stddf = df.groupby(['zbin', 'statname']).std()
    else:
        stddf = df.groupby('statname').std()

    assert stddf.filter(like='bin_').eq(0).all().all(), errmsg
    return None


# look at galaxy distribution
def plot_galaxies(galaxies, gal_frac=0.05, title='Galaxies', coords='xyz', save=None):
    """ Plots 3D galaxy distribution.
        galaxies assumed to be DataFrame with minimum columns {'x','y','z'}.
            Must include 'Redshift' column if coords='rz'.
            Including column 'zbin' will use this to color points.
        gal_frac NOT CURRENTLY USED!
        coords = 'xyz' plots all 3 coordinates (x,y,z)
                 'rz' plots sqrt(x^2+y^2+z^2) vs. Redshift (useful to see redshift binning)
        save = string. path to save plot
    """
    # get a random sample
    l = len(galaxies)
    max_points = 5e4 # max number of points to plot
    gal_frac = min(1,max_points/l) # fraction of points to actually plot
    gs = galaxies.sample(int(l*gal_frac))

    plt.figure()
    proj = '3d' if coords=='xyz' else None
    ax = plt.axes(projection=proj)
    # plt.figure().gca(projection='3d')

    # plot and color by zbin if available
    c = gs.zbin if ('zbin' in gs.columns) else 'b'

    if coords=='xyz':
        ax.scatter3D(gs.x, gs.y, gs.z, s=1, c=c)
        ax.set_xlabel(r'$x\ h^{-1}$ [Mpc]')
        ax.set_ylabel(r'$y\ h^{-1}$ [Mpc]')
        ax.set_zlabel(r'$z\ h^{-1}$ [Mpc]')
    elif coords=='rz':
        # sz = c*1000
        r = np.sqrt(gs.x**2 + gs.y**2 + gs.z**2)
        ax.scatter(r, gs.Redshift, s=1, alpha=0.5, c=c)
        # plt.loglog()
        plt.xlabel(r'r\ h^{-1}$ [Mpc]')
        plt.ylabel('Redshift')

        # calc and display num galaxies in each zbin
        str = 'zbin   # Galaxies'
        for i, (zzz,gsgz) in enumerate(gs.groupby('zbin')):
            str = str+ '\n{z:4.2f}  {gals}'.format(z=zzz, gals=len(gsgz))
        plt.annotate(str, (0.15,0.75), xycoords='figure fraction')

    else: # raise an error
        assert 0, 'plot_galaxies() received invalid argument coords = {}'.format(coords)

    plt.title(title)
    plt.tight_layout()
    if save is not None:
        plt.savefig(save)
    plt.show(block=False)
    plt.pause(10.0)

    return None
