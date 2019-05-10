import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import FormatStrFormatter
# import imp

from astropy import cosmology
import calc_wtheta as cw
import helper_fncs as hf

# set plot defaults
mpl.rcParams['figure.figsize'] = [14.0, 8.0]
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
    validate_statmeans_xbins(df) # make sure we haven't averaged different xbins
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


def get_bins_stats(row, stat, avg_zbins=False):
    """ Returns 2 series plus a y-axis label:
        x = columns starting with 'bin_'
        y = columns starting with 'stat_'
        Scales y values according to stat.

        if avg_zbins==True, row should be a DataFrame
    """

    lbldict = { 'wtheta': r'$\theta\ w(\theta)$', \
                'wp': r'$r\ w_p(r_p)$', \
                'xi': r'$r^2\ \xi(r)$'}

    if not avg_zbins:
        x = row.filter(like='bin_').rename(index=lambda xx: xx.split('_')[-1])
        y = row.filter(like='stat_').rename(index=lambda xx: xx.split('_')[-1])
        if stat in ['wtheta', 'wp']:
            y = x*y
        elif stat in ['xi']:
            y = x*x*y
        else:
            print('stat {} not listed in myplots.get_bins_stats()'.format(stat))

    else:
        assert stat=='wtheta', 'mp.get_bins_stats() received avg_zbins=True but stat!=wtheta'
        # Get theta*wtheta, then take average
        x = row.filter(like='bin_').rename(columns=lambda xx: xx.split('_')[-1]) # df
        y = row.filter(like='stat_').rename(columns=lambda xx: xx.split('_')[-1]) # df
        y = x*y # df
        x = x.mean(axis=0) # series
        y = y.mean(axis=0) # series

    return (x, y, lbldict[stat])


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


def plot_wtheta(fdat, spcols = ['Nstack','NR/NG'], avg_zbins=False, save=None, show=True):
    """ Plots wtheta(theta_bcens)
        fdat (string): path to stats.dat file as written by MockBox.write_stat_to_file()
        spcols (list): ['row_name','col_name']. Does 'groupby' on wdf to plot unique row_name
                    and col_name values in subplots.
        avg_zbins (bool): ==    False plots one line for each zbin
                                True plots one line with wtheta averaged over z
                                        (for use with MB.tratio_binedges)
    """

    wdf = load_statsdat(fdat, stat='wtheta')
    validate_statmeans_xbins(wdf, avg_zbins=avg_zbins)
        # make sure we don't averaged different theta bins

    # Set up subplots
    if 'NR/NG' in spcols: # create this column in wdf
        wdf['NR/NG'] = np.round((wdf['Nrands']/wdf['Ngals']),2)
    rcol, ccol = spcols[0], spcols[1]
    nrows, ncols = len(wdf[rcol].unique()), len(wdf[ccol].unique())
    fig, axs = plt.subplots(nrows, ncols, sharex=True, sharey=True)

    # 2 groupby's to get df for individual subplots
    wdf_rowgroup = wdf.groupby(rcol)
    for i, (rkey, rdf) in enumerate(wdf_rowgroup):
        # wdf_rowcolgroup = rdf.groupby(ccol).mean() # df
        # lendf = rdf.groupby(ccol).size() # series with # of mocks aggragated in each df above
        # for j, (ckey, row) in enumerate(wdf_rowcolgroup.iterrows()):
        wdf_rowcolgroup = rdf.groupby(ccol)
        for j, (ckey, rcdf) in enumerate(wdf_rowcolgroup):
            if (nrows == 1) & (ncols == 1):
                ax = axs
            elif nrows == 1:
                ax = axs[j]
            elif ncols == 1:
                ax = axs[i]
            else:
                ax = axs[i,j]

            # Do groupbys and get the data to plot
            if not avg_zbins:
                rcdfg = rcdf.groupby('zbin').mean() # DF, index = zbin
                lendf = rcdf.groupby('zbin').size() # series with # of mocks aggragated in each df above
                zpm = row.zwidth/2
                for zzz, row in rcdfg.iterrows():
                    x,y,ylabel = get_bins_stats(row, 'wtheta')
                    lbl = r'{:7.2f}$\pm${:.2f} {:5.0f} {:5.0f} {:5.0f}'.format(\
                                    zzz, zpm, row.Nstack, lendf.loc[zzz], row['NR/NG'])
                    ax.semilogx(x,y, label=lbl)
            else: # do same as above but don't groupby zbin
                x,y,ylabel = get_bins_stats(rcdf, 'wtheta', avg_zbins=avg_zbins)
                zzz = rcdf.zbin.mean()
                lendf = len(rcdf)
                zpm = (rcdf.zbin.max()-rcdf.zbin.min()+rcdf.zwidth.mean())/2
                lbl = r'{:7.2f}$\pm${:.2f} {:5.0f} {:5.0f} {:5.0f}'.format(\
                                zzz, zpm, rcdf.Nstack.mean(), lendf, rcdf['NR/NG'].mean())
                ax.semilogx(x,y, label=lbl)

            ax.axhline(0, c='0.5')
            ax.grid(which='both')
            # ax.set_title('wtheta')
            zlbl = '                  zbin' if not avg_zbins else '                  zavg'
            ax.legend(title='{:25s}{:6s}{:6s}{:4s}'.format(zlbl, 'Nstk', '#Moks', 'R/G'), loc=3)
            # ax.set_xlabel(r'$\theta$ [deg]' if not avg_zbins else r'$\theta_{avg}$ [deg]')
            # ax.set_ylabel(ylabel)

            # if stat != 'wtheta':
            #     ax.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
            #     ax.xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))


            # Title subplots with rkey, ckey
            if i==0: # top row
                ax.set_title('{colname} = {colkey}'.format(colname=ccol, colkey=ckey))
            elif i==nrows-1: # bottom row
                ax.set_xlabel(r'$\theta$ [deg]' if not avg_zbins else r'$\theta_{avg}$ [deg]')
            if j==0: # left column
                ax.set_ylabel(ylabel)
            if j==ncols-1: # right column
                ax.annotate('{rowname} = {rowkey}'.format(rowname=rcol, rowkey=rkey), \
                                (1.05,0.5), xycoords='axes fraction', rotation=-90)

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.tight_layout()
    if save is not None:
        plt.savefig(save)
    if show:
        plt.show(block=False)

    return None



def load_statsdat(fdat, stat=None, clean=True):
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




def plot_wtheta_old(wdf, spcols = ['Nstack','NR/NG'], save=None, show=True):
    """wdf = DataFrame with columns wtheta(theta_bcens), 'zbin', 'mock'
        (if multiple mock nums for given zbin, get average.)
        Assumes all column names that can be converted to a float are theta bin centers
        spcols = list ['row_name','col_name']. Does 'groupby' on wdf to plot unique row_name
                    and col_name values in subplots.
        Plots wtheta(theta_bcens), one line for each zbin
    """

    bincols, ocols = cw.get_tbins(wdf) # get list of theta bin and other column names

    # Set up subplots
    if 'NR/NG' in spcols: # create this column in wdf
        wdf['NR/NG'] = (wdf['Nrands']/wdf['Ngals']).astype(int)
    rcol, ccol = spcols[0], spcols[1]
    nrows, ncols = len(wdf[rcol].unique()), len(wdf[ccol].unique())
    fig, axs = plt.subplots(nrows, ncols, sharex=True, sharey=True)
    plt.ylim(-0.002,0.006)
    plt.xlim(0.7,max([float(b) for b in bincols])+1.)

    # 2 groupby's to get df for individual subplots
    wdf_rowgroup = wdf.groupby(rcol)
    for i, (rkey, rdf) in enumerate(wdf_rowgroup):
        # print('i=',i)
        wdf_rowcolgroup = rdf.groupby(ccol)
        for j, (ckey, df) in enumerate(wdf_rowcolgroup):
            # print('j=',j)
            # Create new df of mean wtheta values for each zbin using pivot_table.
                # Transpose to use df.plot()
            wtheta = (pd.pivot_table(df[bincols+['zbin']], index='zbin')).T
                    # cols = zbin, rows = thetabin
            wtheta.rename(index=lambda c: np.double(c), inplace=True)
                    # change index dtype to double

            # Draw the subplot
            if (nrows == 1) & (ncols == 1):
                ax = axs
            elif nrows == 1:
                ax = axs[j]
            elif ncols == 1:
                ax = axs[i]
            else:
                ax = axs[i,j]
            wtheta.sort_index().plot(ax=ax)
            ax.axhline(0, c='0.5')

            # Annotate with extra info
            nmocks = pd.pivot_table(df[['zbin','mock']], index='zbin', aggfunc=len)
                # count number of mocks in the average
            ngals = pd.pivot_table(df[['zbin','Ngals']], index='zbin')
            ngals_std = pd.pivot_table(df[['zbin','Ngals']], index='zbin', aggfunc=np.std)
            ngals_std.rename(columns={'Ngals':'Ngals_std'}, inplace=True)
            str = '{z:^9s} {nm:^10s} {ang:^22s}'.format(z='zbin', nm='Nmocks', ang='Ngals')
            for z in wtheta.columns.values:
                nm = nmocks.loc[z].mock
                ang = ngals.loc[z].Ngals
                ngstd = ngals_std.loc[z].Ngals_std if nm!=1 else 0
                str = str+ '\n{z:^9.2f} {nm:^10.0f} {ang:13.1e}$\pm${ngstd:<9.0e}'.format(\
                        z=z, nm=nm, ang=ang, ngstd=ngstd)
            ax.annotate(str, (0.4,0.75), xycoords='axes fraction')

            # Title subplots with rkey, ckey
            if i==0: # top row
                ax.set_title('{colname} = {colkey}'.format(colname=ccol, colkey=ckey))
            elif i==len(wdf_rowgroup.groups)-1: # bottom row
                ax.set_xlabel(r'$\theta$ [deg]')
            if j==0: # left column
                ax.set_ylabel(r'$w(\theta)$')
            else: # right column (assumes 2 columns)
                # ax.set_ylabel('{rowname} = {rowkey}'.format(rowname=rcol, rowkey=rkey))
                ax.annotate('{rowname} = {rowkey}'.format(rowname=rcol, rowkey=rkey), \
                                (-0.05,0.9), xycoords='axes fraction')
            if ncols==1:
                ax.annotate('{rowname} = {rowkey}'.format(rowname=rcol, rowkey=rkey), \
                                (1,0.9), xycoords='axes fraction', rotation=-90)


    # plt.title('Average of {:.1f} mocks'.format(len(wdf)/len(wdf.zbin.unique())))
    plt.semilogx()
    ax.xaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
    plt.tight_layout()
    if save is not None:
        # plt.ylim(-0.0015, 0.002)
        # plt.tight_layout()
        plt.savefig(save)
    if show:
        plt.show(block=False)

    return None

def plot_wtheta_old_old(wdf, save=None, show=True):
    """wdf = DataFrame with columns wtheta(theta_bcens), 'zbin', 'mock'
        (if multiple mock nums for given zbin, get average.)
        Assumes all column names that can be converted to a float are theta bin centers
        Plots wtheta(theta_bcens), one line for each zbin
    """

    bincols, ocols = cw.get_tbins(wdf) # get list of theta bin column names
    # create new df of wtheta values using pivot_table
        # calculating mean wtheta for each zbin.
        # Transpose to use df.plot()
    wdf.zbin = np.round(wdf.zbin,2) # should store zbin this way from the beginning
    wtheta = (pd.pivot_table(wdf[bincols+['zbin']], index='zbin')).T # cols = zbin, rows = thetabin
    wtheta.rename(index=lambda c: np.double(c), inplace=True) # change index dtype to double
    wtheta = wtheta.sort_index()

    plt.figure()
    # plt.scatter()
    wtheta.plot()
    # wtheta.plot(x = , y = , kind='scatter')
    plt.axhline(0, c='0.5')
    # plt.scatter(wtheta.index.values, wtheta.)

    # annotate with other info for each zbin
    str = 'Nstack avg = {:2.1f}\nzbin  Ngals avg  NR/NG avg'.format(np.average(wdf.Nstack.unique()))
    desc_df = pd.pivot_table(wdf[ocols], index='zbin')
    ddfg = desc_df.groupby('zbin')
    for i, (zzz, ddf) in enumerate(ddfg):
        str = str+ '\n{:4.2f} {:11.2e} {:11.2f}'.format(zzz, \
                np.average(ddf.Ngals.values), np.average(ddf.Nrands.values/ddf.Ngals.values))
    plt.annotate(str, (0.4,0.75), xycoords='axes fraction')

    plt.xlabel(r'$\theta$ [deg]')
    plt.ylabel(r'$w(\theta)$')
    plt.title('Average of {:.1f} mocks'.format(len(wdf)/len(wdf.zbin.unique())))
    if save is not None:
        plt.xlim(-0.0015, 0.002)
        plt.tight_layout()
        plt.savefig(save)
    if show:
        plt.show(block=False)

def plot_wtheta_old_old_old(bcens, wtheta):
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
    max_points = 1e5 # max number of points to plot
    gal_frac = min(1,max_points/l) # fraction of points to actually plot
    gs = galaxies.sample(int(l*gal_frac))
    # lg = np.arange(len(galaxies['x']))
    # np.random.shuffle(lg)
    # lg = lg[:int(gal_frac*len(galaxies['x']))]

    plt.figure()
    proj = '3d' if coords=='xyz' else None
    ax = plt.axes(projection=proj)
    # plt.figure().gca(projection='3d')

    # plot and color by zbin if available
    c = gs.zbin if ('zbin' in gs.columns) else 'b'

    if coords=='xyz':
        ax.scatter3D(gs.x, gs.y, gs.z, s=1, c=c)
        ax.set_xlabel(r'x [h$^{-1}$ Mpc]')
        ax.set_ylabel(r'y [h$^{-1}$ Mpc]')
        ax.set_zlabel(r'z [h$^{-1}$ Mpc]')
    elif coords=='rz':
        # sz = c*1000
        r = np.sqrt(gs.x**2 + gs.y**2 + gs.z**2)
        ax.scatter(r, gs.Redshift, s=1, alpha=0.5, c=c)
        # plt.loglog()
        plt.xlabel(r'r [h$^{-1}$ Mpc]')
        plt.ylabel('Redshift')

        # calc and display num galaxies in each zbin
        gsg = gs.groupby('zbin')
        str = 'zbin   # Galaxies'
        for i, (zzz,gsgz) in enumerate(gsg):
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


def plot_dtheta_rp(zlist=None, tratio_binedges=None, save=None):
    """ Plots delta theta [deg] vs projected distance [h^-1 Mpc]
        for different redshifts.

        Calculates tbin_edges using tratio_binedges for each z in zlist.
        Then calculates rp (projected distance) for each theta.
    """

    if zlist is None:
        zlist = [0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]

    if tratio_binedges is None:
        tratio_binedges = np.logspace(np.log10(0.1), np.log10(2.), 11)

    plt.figure()
    cmap = plt.get_cmap('inferno_r')
    colors = [cmap(c) for c in np.linspace(0.1,0.9,len(zlist)+1)]
    for i, z in enumerate(zlist):
        t_binedges, rp_binedges = hf.get_theta_rp_from_tratio_bins(z, tratio_binedges)
        lbl = '{}'.format(z)
        plt.loglog(t_binedges, rp_binedges, '+-', ms=6, c=colors[i], label=lbl)
        if np.mod(i,2) == 0:
            plt.plot(t_binedges[-1], rp_binedges[-1], 'k+')
            plt.annotate(lbl, (t_binedges[-1], rp_binedges[-1]+10))

    plt.axhline(105, c='k')
    plt.legend(title='Redshift')
    plt.xlabel(r'$\Delta \theta$ [deg]')
    plt.ylabel(r'projected distance $h^{-1}$[Mpc]')
    plt.grid(axis='x', which='both')
    if save is not None:
        plt.savefig(save)
    plt.show(block=False)

    return None




def plot_redshift_thetaBAO(cosmo=None, save=None):
    """ plots expected theta of BAO feature (mapped from d_BAO = 150 Mpc)
            as a function of redshift.
        save = string of path to save plot
    """

    if cosmo is None:
        cosmo = cosmology.FlatLambdaCDM(H0=70.0, Om0=0.3)

    zz = np.logspace(np.log10(0.1), np.log10(3.0), int(1e4))
    rz = cosmo.comoving_distance(zz).value # = r(zz) [Mpc]
    rzh = cosmo.comoving_distance(zz).value*cosmo.h # = r(zz) [Mpc/h]
    dBAO = 150 # [Mpc/h]
    thetaBAO = np.degrees(np.arctan2(dBAO,rz))
    thetaBAOh = np.degrees(np.arctan2(dBAO,rzh))

    ax = plt.subplot(111)
    plt.plot(zz, thetaBAO, label='r(z) in [Mpc]')
    plt.plot(zz, thetaBAOh, label='r(z) in [Mpc/h]')
    plt.semilogy()
    plt.xlabel('Redshift')
    plt.ylabel(r'$\theta_{BAO} = atan(\frac{150Mpc}{r(Redshift)})$')
    plt.grid(which='both')
    plt.ylim(1.2,13)
    ax.yaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.0f"))
    plt.legend()
    plt.title(r'Expected $\theta$ of BAO signal in w($\theta$)')
    if save is not None:
        plt.savefig(save)
    plt.show(block=False)

    return None


# plot comoving distance vs redshift
# call with:
    # halocat, HODmodel = sm.setup_halosHOD() # fetch halo catalog and HODmodel,
        # returns populated HODmodel.mock
    # zred = halocat.redshift
    # plot_codist_redshift(log=False, fout=False, zred=zred)
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
