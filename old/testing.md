# Use Outer Rim Fof halos
get metadata: Lbox [Mpc/h], particle_mass [Msun/h], redshift, simname='outerrim'
get halo data (arrays): halo_id, halo_x, halo_y, halo_z, halo_mass
halo_catalog = UserSuppliedHaloCatalog(redshift = redshift, Lbox = Lbox, particle_mass = particle_mass, halo_x = halo_x, halo_y = halo_y, halo_z = halo_z, halo_id = halo_id, halo_mass = halo_mass)




# Fit correlation functions
<!-- fs -->
```python
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import myplots as mp
from scipy.optimize import curve_fit
def corr_fitfnc(x, A, gamma, x0):
    return A* (x/x0)**(-gamma)
fdat = 'data/stats_tratiobins_zw0.15.dat'
df = mp.load_statsdat(fdat)
xdf = df.query("statname=='{}'".format('xi')).mean()
# x,y,lbl = mp.get_bins_stats(xdf, 'xi')
x = xdf.filter(like='bin_')
x.index = np.arange(len(x))
y = xdf.filter(like='stat_')
y.index = np.arange(len(y))
plt.figure()
plt.scatter(x,x*x*y)
plt.plot(x, x*x*corr_fitfnc(x, 1., 1.8, 5.0))
plt.show(block=False)
```
<!-- fe # Fit correlation functions -->


# Test 2: tratio_binedges
<!-- fs -->
Using wider theta bins avoids NaNs (see sections [here](#deal-with-nans-in-wtheta-stats-output) and [here](#test-1-tratio_binedges))

`python -u main.py >> maintratiobins2.out`
<!-- fs run main.py with:
# DEFAULTS:
statfout='data/stats.dat'
stats=['wtheta', 'xi', 'wp']
nbins = 50
# tbin_edges = np.logspace(np.log10(0.05), np.log10(9.0), nbins+1)
tratio_binedges = np.logspace(np.log10(0.05), np.log10(2.), nbins+1)
rbin_edges = np.logspace(np.log10(25.0), np.log10(150.0), nbins+1)
pimax = 300
galplots = False
z4push = 'cat'
zw_list = [0.1, 0.05, 0.15]
cat_gals = 5e5 # approx num gals in cat mock
Nstack_list = [0, 2] # number of mock boxes to stack, per dimension
nrfact_list = [1, 3, 10] # used in Nrands
# Nrands = int(nrfact*cat_gals* max(Nstack,1)**3)
imax = 5 # number of times to run each param combo
 -->
<!-- Plots:
import myplots as mp
tags = ['stats_tratiobins_zw0.1', 'stats_tratiobins_zw0.05', 'stats_tratiobins_zw0.15']
for tag in tags:
    mp.plot_wtheta('data/'+tag+'.dat', avg_zbins=True, save='plots/'+tag+'.png')
 -->
<!-- fe run main.py with: -->

- [x] zbin width = 0.05.
    - tag: stats_tratiobins_zw0.05
    <img src="plots/stats_tratiobins_zw0.05.png" alt="stats_tratiobins_zw0.05" width="800"/>

- [x] zbin width = 0.1.
    - tag: stats_tratiobins_zw0.1
    <img src="plots/stats_tratiobins_zw0.1.png" alt="stats_tratiobins_zw0.1" width="800"/>

- [x] zbin width = 0.15.
    - tag: stats_tratiobins_zw0.15
    <img src="plots/stats_tratiobins_zw0.15.png" alt="stats_tratiobins_zw0.15" width="800"/>


<!-- fe # Test tratio_binedges -->


# Deal with NaNs in wtheta stats output
<!-- fs -->
Corrfunc returns NaN from wtheta calculation for "bins where the RR count is 0" (documentation page 85).

```python
import numpy as np
import myplots as mp
import matplotlib as mpl
from matplotlib import pyplot as plt
fsog = mpl.rcParams['figure.figsize']
mpl.rcParams['figure.figsize'] = [14.0, 3.0]
# mpl.rcParams['figure.figsize'] = fsog

taglst = ['stats_tratiobins_zw0.05', 'stats_tratiobins_zw0.1', 'stats_tratiobins_zw0.15']
flist = ['data/'+tag+'.dat' for tag in taglst]
df = mp.load_statsdat(flist, stat='wtheta', clean=True)
df['RG'] = np.round(df['Nrands']/df['Ngals'], 2) # NR/NG

# Plot total number of NaNs per column
nuls = df.isnull().sum(axis=0)
nuls[nuls>0].plot()
plt.ylabel('# NaN values');
plt.title('stats_tratiobins_zw[0.05, 0.1, 0.15]; {} total rows'.format(len(df)));
plt.savefig('plots/test/stat_Nans_per_thetabin.png'); plt.show(block=False)

# # cleaup for groupbys
# df['zwidth'] = 0.05*np.round(df.zwidth/0.05)
# df.loc[df.zwidth>0.15,'zwidth'] = 0.15
# df.loc[df.RG==4, 'RG'] = 3; df.loc[df.RG>10, 'RG'] = 10
# # Plot total number of NaNs per {zwidth, Nstack, RG}
# sz = df.groupby([df.Nstack, df.zwidth, df.RG]).size()
# nas = df.isnull().sum(axis=1).groupby([df.Nstack, df.zwidth, df.RG]).sum()/sz
# nas.plot(marker='o');
# plt.xticks(np.arange(len(nas)), rotation=60); plt.gca().set_xticklabels(nas.index.values);
# plt.ylabel('# NaN values'); plt.tight_layout(); plt.show(block=False)

# Plot average number of NaNs per zbin
df['zwidth'] = 0.05*np.round(df.zwidth/0.05)
dg = df.groupby('zwidth')
plt.figure()
for zw, d in dg:
    sz = d.groupby('zbin').size()
    nas = d.isnull().sum(axis=1).groupby(d.zbin).sum()/sz
    nas.plot(marker='o', label='zwidth = {}'.format(zw));
plt.legend(); plt.ylabel('avg # NaNs'); plt.tight_layout();
plt.savefig('plots/test/stat_avgNaNs_per_zbin.png'); plt.show(block=False)

# Plot average number of NaNs per RG
dg = df.groupby('Nstack')
plt.figure()
for nst, d in dg:
    sz = d.groupby('RG').size()
    nas = d.isnull().sum(axis=1).groupby(d.RG).sum()/sz
    plt.scatter(nas.index.values, nas, label='Nstack = {}'.format(nst))
plt.legend(); plt.ylabel('avg # NaNs'); plt.xlabel('# randoms per galaxy in sample');
plt.tight_layout(); plt.show(block=False)
plt.savefig('plots/test/stat_avgNaNs_per_NRNG.png');
# nacols = df.isnull().any(axis=0)
# narows = df.isnull().any(axis=1)
```
- [x] Drop rows of wdf where Nrands < 1000
    - check zbin of these rows. Done: they are all the last (max) zbin as expected

- [x] Total \# NaN values in each theta bin ('stat_#' column)
    - there are 75 theta bins, all those not shown have 0 NaNs
    - all NaNs are in the 15 smallest theta bins => should try starting bins at larger theta and making each bin larger.
    <img src="plots/test/stat1_Nans_per_thetabin.png" alt="stat1_Nans_per_thetabin" width="800"/>

- [x] Average \# NaN values in each redshift bin
    - Theta bins correspond to constant projected distance, independent of redshift
    - As redshift increases, each zbin width corresponds to a smaller comoving distance(width)
    <img src="plots/test/stat1_avgNaNs_per_zbin.png" alt="stat1_avgNaNs_per_zbin" width="800"/>

- [x] Average \# NaN values as fnc of \# randoms per galaxy in sample
    - \# NaNs should decrease as NR/NG increases
    <img src="plots/test/stat1_avgNaNs_per_NRNG.png" alt="stat1_avgNaNs_per_NRNG" width="800"/>


<!-- fe # Deal with NaNs in wtheta stats output -->


# Test 1: tratio_binedges
<!-- fs -->
`python -u main.py >> maintratiobins.out`
<!-- fs run main.py with:

# DEFAULTS:
stats=['wtheta', 'xi', 'wp']
nbins = 75
# tbin_edges = np.logspace(np.log10(0.05), np.log10(9.0), nbins+1)
tratio_binedges = np.logspace(np.log10(0.01), np.log10(2.), nbins+1)
rbin_edges = np.logspace(np.log10(25.0), np.log10(150.0), nbins+1)
pimax = 300
galplots = False
z4push = 'cat'
zw_list = [0.1, 0.05, 0.15]
cat_gals = 5e5 # approx num gals in cat mock
Nstack_list = [0, 2] # number of mock boxes to stack, per dimension
nrfact_list = [1, 3, 10] # used in Nrands
# Nrands = int(nrfact*cat_gals* max(Nstack,1)**3)
imax = 5 # number of times to run each param combo
 -->
<!-- Plots:
import myplots as mp
tags = ['stats_tratiobins_zw0.1', 'stats_tratiobins_zw0.05', 'stats_tratiobins_zw0.15']
for tag in tags:
    mp.plot_wtheta('data/'+tag+'.dat', avg_zbins=True, save='plots/'+tag+'.png')
 -->
<!-- fe run main.py with: -->

- [x] zbin width = 0.05. (ow_051019_1557)
    - tag: stats_tratiobins_zw0.05
    <img src="plots/stats_tratiobins1_zw0.05.png" alt="stats_tratiobins1_zw0.05" width="800"/>

- [x] zbin width = 0.1.
    - tag: stats_tratiobins_zw0.1 (ow_051019_1341)
    <img src="plots/stats_tratiobins1_zw0.1.png" alt="stats_tratiobins1_zw0.1" width="800"/>

- [x] zbin width = 0.15. (ow_051019_1742)
    - tag: stats_tratiobins_zw0.15
    <img src="plots/stats_tratiobins1_zw0.15.png" alt="stats_tratiobins1_zw0.15" width="800"/>


<!-- fe # Test tratio_binedges -->


# [x] Test NR/NG = 0.1, 0.01
<!-- fs -->
`python -u main.py >> mainNRNGsmall.out`
<!-- fs run main.py with:

# DEFAULTS:
statfout='data/stats.dat'
stats=['wtheta', 'xi', 'wp']
nbins = 75
tbin_edges = np.logspace(np.log10(0.05), np.log10(9.0), nbins+1)
rbin_edges = np.logspace(np.log10(25.0), np.log10(150.0), nbins+1)
pimax = 300
galplots = False
z4push = 'cat'
zw_list = [10, 0.3]
cat_gals = 5e5 # approx num gals in cat mock
Nstack_list = [0, 2] # number of mock boxes to stack, per dimension
nrfact_list = [0.01, 0.1] # used in Nrands
# Nrands = int(nrfact*cat_gals* max(Nstack,1)**3)
imax = 20 # number of times to run each param combo
 -->
<!-- Plots:
import myplots as mp
statfout='data/stats_NRNGsmall_zw10.dat'
mp.plot_wtheta(statfout, save='plots/stats_NRNGsmall_zw10.png')
statfout='data/stats_NRNGsmall_zw0.3.dat'
mp.plot_wtheta(statfout, save='plots/stats_NRNGsmall_zw0.3.png')
 -->

<!-- fe fs python -u main.py >> mainNRNGsmall.out -->
- [x] Full box (no redshift binning). tag: stats_NRNGsmall_zw10
    <img src="plots/stats_NRNGsmall_zw10.png" alt="stats_NRNGsmall_zw10" width="800"/>

- [x] Redshift bin width = 0.3. tag: stats_NRNGsmall_zw0.3
    <img src="plots/stats_NRNGsmall_zw0.3.png" alt="stats_NRNGsmall_zw0.3" width="800"/>

<!-- fe # Test NR/NG = 0.1, 0.01 -->



# [x] Test Corrfunc.theory xi, wp
<!-- fs -->
<!-- fs Run main.py with:
python -c "import helper_fncs as hf; hf.file_ow('main.out')"
python -u main.py >> main.out # -u forces unbuffered stdout

# DEFAULTS:
statfout='data/stats.dat'
stats=['wtheta', 'xi', 'wp']
nbins = 51
tbin_edges = np.logspace(np.log10(0.1), np.log10(12.0), nbins+1)
rbin_edges = np.logspace(np.log10(75.0), np.log10(150.0), nbins+1)
pimax = 300
galplots=False
z4push = 'cat'
zw = 10.
cat_gals = 5e5 # approx num gals in cat mock
Nstack = 2 # number of mock boxes to stack, per dimension
nrfact = 1 # used in Nrands
Nrands = int(nrfact*cat_gals* max(Nstack,1)**3)
imax = 10 # number of times to run each param combo
-->
<!-- Plots:
import myplots as mp
statfout='data/stats_Nstack0_z4push0.dat'
mp.plot_stats(statfout, save='plots/stats_Nstack0_z4push0.png', show=True)
statfout='data/stats_Nstack0.dat'
mp.plot_stats(statfout, save='plots/stats_Nstack0.png', show=True)
statfout='data/stats_defaults.dat'
mp.plot_stats(statfout, save='plots/stats_defaults.png', show=True)
-->

<!-- fe Run main.py with: -->

- [x] Nstack=0, z4push=0
    - statfout='data/stats_Nstack0_z4push0.dat'
    <img src="plots/stats_Nstack0_z4push0.png" alt="stats_Nstack0_z4push0" width="800"/>


- [x] Nstack=0
    - statfout='data/stats_Nstack0.dat'
    <img src="plots/stats_Nstack0.png" alt="stats_Nstack0" width="800"/>


- [x] Defaults
    - statfout='data/stats_defaults.dat'
    <img src="plots/stats_defaults.png" alt="stats_defaults" width="800"/>



<!-- fe # Test Corrfunc.theory xi, wp -->
