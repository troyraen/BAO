from warnings import warn as _warn
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

cosmo = None

def plot_stats(fdat, param_dict, save=None, show=True, zbin='avg'):
    """ Plots stats (1 per column) in file fdat.
        Args:
        fdat (string): path to stats.dat file as written by MockBox.write_stat_to_file()
        param_dict   : only need key = 'cosmo', value = astropy cosmology object.
                       Needed for wtheta x-axis conversion.
        zbin         : == 'avg': averages zbins in wtheta plot
                       == n (int [0, num zbins - 1]) plots wtheta for nth zbin only
    """
    globals()['cosmo'] = param_dict['cosmo'] # needed for wtheta x-axis conversion

    holdfsz = mpl.rcParams['figure.figsize'] # keep to reset later
    mpl.rcParams['figure.figsize'] = [14.0, 4.0]

    df = load_statsdat(fdat)
    if type(zbin) == int: # keep nth zbin
        zbin_cens = np.sort(df.loc[df['statname']=='wtheta','zbin'].unique())
        zzz = zbin_cens[zbin]
        df.drop(labels=df.loc[((df['statname']=='wtheta') & (df['zbin']!=zzz))].index,
                inplace=True)

    sdf = df.groupby('statname').mean() # df
    sdf['NR/NG'] = (sdf['Nrands']/sdf['Ngals']).astype(int)
    # validate_statmeans_xbins(df) # make sure we haven't averaged different xbins
    lendf = df.groupby('statname').size() # series with # of mocks aggragated in each df above

    nrows, ncols = 1, len(lendf)
    fig, axs = plt.subplots(nrows, ncols, sharex=False, sharey=False)

    # Plot for each stat
    for i, (stat, row) in enumerate(sdf.iterrows()):
        ax = axs[i]

        if zbin == 'avg' and stat == 'wtheta':
            x,y,axlbls = get_bins_stats(df.loc[df.statname=='wtheta',:], stat, avg_zbins=True)
        else:
            x,y,axlbls = get_bins_stats(row, stat)
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
        ax.set_xlabel(axlbls[0])
        ax.set_ylabel(axlbls[1])

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


def get_bins_stats(row, stat, avg_zbins=False):
    """ Returns 2 series plus a x- and y-axis labels.
        x = columns starting with 'bin_'
        y = columns starting with 'stat_'
        Scales y values according to stat.

        if avg_zbins==True, row should be a DataFrame
    """

    # {statname: (xlabel, ylabel)}
    lbldict = { 'wtheta': (r'$\theta/\theta_{BAO}$', r'$\theta\ w(\theta)$'),
                'wp':     (r'$r_p\ h^{-1}$ [Mpc]', r'$r_p\ w_p(r_p)$'),
                'xi':     (r'$r\ h^{-1}$ [Mpc]', r'$r^2\ \xi(r)$')
              }

    if not avg_zbins:
        x = row.filter(like='bin_').rename(index=lambda xx: xx.split('_')[-1])
        y = row.filter(like='stat_').rename(index=lambda xx: xx.split('_')[-1])
    else:
        x = row.filter(like='bin_').rename(columns=lambda xx: xx.split('_')[-1])
        y = row.filter(like='stat_').rename(columns=lambda xx: xx.split('_')[-1])

    if stat == 'wtheta':
        y = x*y
        if avg_zbins:
            y = y.mean(axis=0)
            x = get_theta_rp_from_tratio_bins(row.zbin.mean(), x.mean(axis=0), invert=True)
        else:
            x = get_theta_rp_from_tratio_bins(row.zbin, x, invert=True)

    elif stat == 'wp':
        y = x*y

    elif stat == 'xi':
        y = x*x*y

    else:
        _warn(f'stat {stat} not listed in myplots.get_bins_stats()')

    return (x, y, lbldict[stat])


# copied from calc_stats.py
def get_theta_rp_from_tratio_bins(redshift, tratio_binedges, invert=False):
    """ Converts tratio_binedges from units of theta_BAO(zbin) to degrees
        at given redshift using law of cosines.

        invert == True: 2nd arg (tratio_binedges) should be t_binedges
                        Returns tratio_binedges
    """

    d_BAO = 105 # h^-1 Mpc
    rz = (cosmo.comoving_distance(redshift).value)*cosmo.h # Mpc/h
    theta_BAO = np.degrees(np.arccos(1 - d_BAO**2/2/rz**2)) # law of cosines

    if not invert:
        t_binedges = theta_BAO* tratio_binedges # degrees
        rp_binedges = rz* np.sqrt(2*(1-np.cos(np.radians(t_binedges)))) # projected, comoving distance
        return t_binedges, rp_binedges
    else:
        trat_binedges = tratio_binedges/theta_BAO
        return trat_binedges

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
    plt.set_cmap('tab20')
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
        plt.xlabel(r'$r\ h^{-1}$ [Mpc]')
        plt.ylabel('Redshift')

        # calc and display num galaxies in each zbin
        str = 'zbin   # Galaxies'
        for i, (zzz,gsgz) in enumerate(gs.groupby('zbin')):
            str = str+ '\n{z:4.2f}  {gals}'.format(z=zzz, gals=len(gsgz))
        plt.annotate(str, (0.15,0.15), xycoords='figure fraction')

    else: # raise an error
        assert 0, 'plot_galaxies() received invalid argument coords = {}'.format(coords)

    plt.title(title)
    plt.tight_layout()
    if save is not None:
        plt.savefig(save)
    plt.show(block=False)
    plt.pause(10.0)

    return None
