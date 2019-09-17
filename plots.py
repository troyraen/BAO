from warnings import warn as _warn
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

cosmo = None
d_BAO = 105 # h^-1 Mpc

def plot_stats(fdat, param_dict, save=None, show=True, zbin='avg', keep_zbin=None):
    """ Plots stats (1 per column) in file fdat.
        Args:
        fdat (string): path to stats.dat file as written by MockBox.write_stat_to_file()
        param_dict   : only need key = 'cosmo', value = astropy cosmology object.
                       Needed for wtheta x-axis conversion.
        zbin         : == 'avg': averages zbins in wtheta plot
                       == 'sep': plots zbins separately in wtheta plot
        keep_zbin    : == None: keeps all zbins (with enough Ngals) for wtheta plot
                       == n (int [0, num zbins - 1]) keeps only the nth zbin(s)
    """
    globals()['cosmo'] = param_dict['cosmo'] # needed for wtheta x-axis conversion

    df = load_statsdat(fdat, clean=True)
    plot_dict = get_stats_plot_data(df, zbin, keep_zbin)

    nrows, ncols = 1, len(df.statname.unique())
    fig, axs = plt.subplots(nrows, ncols, sharex=False, sharey=False, figsize=(14.0,4.0))
    plt.set_cmap('YlOrRd')
    plt.clim(vmin=0, vmax=4)

    for i, (stat, data) in enumerate(plot_dict.items()):
        ax = axs[i]

        for j in range(len(data['x'])):
            lbl='lbl'
            ax.semilogx(data['x'][j], data['y'][j], label=lbl)

        ax.axhline(0, c='0.5')
        ax.grid(which='both')
        ax.set_title(stat)
        ax.legend(title='zbin, Nstk, Nm, NR')
        ax.set_xlabel(data['axlbls'][0])
        ax.set_ylabel(data['axlbls'][1])

    plt.tight_layout()
    if save is not None:
        plt.savefig(save)
    if show:
        plt.show(block=False)

    return None


def get_stats_plot_data(df, zbin, keep_zbin):
    """ Returns dict of stat info for plotting.

    Args:
        df          : stats DataFrame as returned by load_statsdat(fdat, clean=True)
        zbin        : as passed to plot_stats() (affects wtheta plot)
        keep_zbin   : as passed to plot_stats() (affects wtheta plot)

    Returns
        plot_dict   : { <statname>:
                        { 'x':          series of x data (*)
                          'y':          series of y data (*)
                          'axlbls':     tuple of x,y axis labels
                          'numMocks':   # of mocks averaged into data (*)
                          'mean_z':     mean redshift of data (*)
                          'z_width':    mean zbin width of data (*)
                        }
                      }
                                        (*): if zbin=='sep', these are lists of same
    """

    wtdf = df.loc[df['statname']=='wtheta',:] # df with just wtheta
    odf = df.loc[df['statname']!='wtheta',:] # df with everything but wtheta
    plot_dict = {}

    # Get odf plot data
    validate_statmeans_xbins(odf)
    odf_means = odf.groupby('statname').mean() # df
    odf_numMocks = odf.groupby('statname').size() # series
    for (stat, row) in odf_means.iterrows():
        x, y, axlbls = get_bins_stats(row, stat)
        plot_dict[stat] = { 'x': [x],
                            'y': [y],
                            'axlbls': axlbls,
                            'numMocks': [odf_numMocks[stat]],
                            'mean_z': [row.zbin],
                            'z_width': [row.zwidth]
                          }

    # Get wtheta plot data
    stat = 'wtheta'
    validate_statmeans_xbins(wtdf, sep_by_zbin=True)

    if keep_zbin is not None: # keep only nth zbin
        keepz = np.sort(wtdf.zbin.unique())[keep_zbin]
        wtdf.drop(labels=wtdf.loc[wtdf['zbin']!=keepz].index, inplace=True)

    if zbin == 'avg':
        x, y, axlbls = get_bins_stats(wtdf, stat, avg_zbins=True)
        x = [x]
        y = [y]
        numMocks = [len(wtdf)]
        mean_z = [wtdf.zbin.mean()]
        z_width = [wtdf.zwidth.mean()]

    elif zbin == 'sep':
        zdf_means = wtdf.groupby('zbin').mean() # df
        zdf_numMocks = wtdf.groupby('zbin').size() # series
        x,y,numMocks,mean_z,z_width = ([] for i in range(5))
        for (z, row) in zdf_means.iterrows():
            xx, yy, axlbls = get_bins_stats(row, stat)
            x.append(xx)
            y.append(yy)
            numMocks.append(zdf_numMocks[z])
            mean_z.append(z)
            z_width.append(row.zwidth)

    else:
        _warn(f'Invalid value for zbin argument')

    plot_dict['wtheta'] = { 'x': x,
                            'y': y,
                            'axlbls': axlbls,
                            'numMocks': numMocks,
                            'mean_z': mean_z,
                            'z_width': z_width
                          }

    return plot_dict


def load_statsdat(fdat, stat=None, clean=False):
    """ Load correlation stats data from file fdat, as written by MockBox.write_stat_to_file.
        Returns DataFrame of file data.

        fdat (string or list of strings):   stats output file path(s)
        stat (string):  value in column statname in fdat
        clean (bool):   == True will drop rows with Ngals<=10000
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
        df = df.query("Ngals>10000")
    return df


def validate_statmeans_xbins(df, sep_by_zbin=False):
    """ Checks that xbins of collapsed statname have std_dev == 0
            and therefore are all the same.
        sep_by_zbin == True allows for each zbin to have unique theta bins
    """

    errmsg = "Operation groupby 'statname' has averaged values in different theta or r bins."
    if sep_by_zbin:
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



def plot_stats_old(fdat, param_dict, save=None, show=True, zbin='avg'):
    """ Plots stats (1 per column) in file fdat.
        Args:
        fdat (string): path to stats.dat file as written by MockBox.write_stat_to_file()
        param_dict   : only need key = 'cosmo', value = astropy cosmology object.
                       Needed for wtheta x-axis conversion.
        zbin         : == 'avg': averages zbins in wtheta plot
                       == 'sep': plots zbins separately in wtheta plot
                       == n (int [0, num zbins - 1]) plots wtheta for nth zbin only
    """
    globals()['cosmo'] = param_dict['cosmo'] # needed for wtheta x-axis conversion

    holdfsz = mpl.rcParams['figure.figsize'] # keep to reset later
    mpl.rcParams['figure.figsize'] = [14.0, 4.0]

    df = load_statsdat(fdat, clean=True)
    zbin_cens = np.sort(df.loc[df['statname']=='wtheta','zbin'].unique())

    if type(zbin) == int: # keep only nth zbin for wtheta
        zzz = zbin_cens[zbin]
        df.drop(labels=df.loc[((df['statname']=='wtheta') & (df['zbin']!=zzz))].index,
                inplace=True)

    dfw = df.loc[df.statname=='wtheta',:]

    sdf = df.groupby('statname').mean() # df
    sdf['NR/NG'] = (sdf['Nrands']/sdf['Ngals']).astype(int)
    # validate_statmeans_xbins(df) # make sure we haven't averaged different xbins
    lendf = df.groupby('statname').size() # series with # of mocks aggragated in each df above

    nrows, ncols = 1, len(lendf)
    fig, axs = plt.subplots(nrows, ncols, sharex=False, sharey=False, figsize=(14.0,4.0))
    plt.set_cmap('YlOrRd')
    plt.clim(vmin=0, vmax=4)

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
