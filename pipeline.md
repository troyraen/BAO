# Example Usage

# SETUP
From terminal:
```bash
cd Documents/BAO
htenv
```

----------------------------------------------------------------
# Run wtheta calculation beginning to end
<!-- fs -->
In ipython:
```python
%run main.py

# --- OR (this is currently in main.py, but in case it changes...)
from MockBox import MockBox as MB
mb = MB()
mb.getmock_calcwtheta(galplots=True)

# --- OR run this step-by-step:
import setup as su
import calc_wtheta as cw
import myplots as mp
import helper_fncs as hf
from MockBox import MockBox as MB

if su.cosmo is None:
    su.load_cosmo() # loads default global cosmo object plus H0, Om0

mb = MB()
mb.cat_galtbl, mb.cat_Lbox, mb.cat_zbox = su.get_galtbl(getas='DF')
mb.transform_mock(box='PhaseSpace') # Sets mb.PhaseSpace
mb.RDZ = hf.get_ra_dec_z(mb.PhaseSpace, usevel=True)
mb.bin_redshifs(box='RDZ', validate=False)
mb.get_randoms()
rgroups = mb.Randoms.groupby('zbin', axis=0)

zgroups = mb.RDZ[['RA','DEC','zbin']].groupby('zbin', axis=0)
for i, (zzz, rdz_z) in enumerate(zgroups):
    Randoms_z = rgroups.get_group(zzz) # get the randoms within this zbin
    mb.calc_write_wtheta(zzz, rdz_z, Randoms_z, fout, nthreads=nthreads)

```

OR in Bash:
```bash
# first move main.out if it exists
python -c 'import helper_fncs as hf; hf.file_ow('main.out')'
python -u main.py >> main.out # -u forces unbuffered stdout
```

<!-- fe wtheta calculation -->


----------------------------------------------------------------
# get and plot wtheta df from file
<!-- fs -->
```python
import pandas as pd
import myplots as mp
import calc_wtheta as cw
fin = 'data/wtheta.dat'
wdf = cw.load_from_file(fin)
fout = 'plots/wtheta.png'
mp.plot_wtheta(wdf, save=fout)

# the rest is not working:
wdfp = pd.pivot_table(wdf, index='zbin')
wdfp
mp.plot_wtheta(wdf)

wdfg = wdf.groupby('zbin')
for i, (zzz, df) in enumerate(wdfg):
    mp.plot_wtheta(df)
```
<!-- fe get, plot wtheta from file -->

----------------------------------------------------------------
# get and plot a galaxy table
<!-- fs -->
```python
# get the table
import helper_fncs as hf
import setup as su
if su.cosmo is None:
    su.load_cosmo() # loads default global cosmo object plus H0, Om0
from MockBox import MockBox as MB
zw = 0.4
mb = MB(Nstack=2, zbin_width=zw)
mb.cat_galtbl, mb.cat_Lbox, mb.cat_zbox = su.get_galtbl(getas='DF')
mb.transform_mock(box='PhaseSpace') # Sets mb.PhaseSpace
mb.RDZ = hf.get_ra_dec_z(mb.PhaseSpace, usevel=True)
mb.bin_redshifs(box='RDZ', validate=False)

# plot the table
import pandas as pd
import myplots as mp
df = pd.concat([mb.PhaseSpace[['x','y','z']], mb.RDZ[['Redshift','zbin']]],axis=1)
# plot r vs redshift
mp.plot_galaxies(df, coords='rz', title='zbin_width = {}'.format(zw), save='plots/gals_rvsz.png')
# plot x,y,z coordinates
mp.plot_galaxies(df, coords='xyz', title='zbin_width = {}'.format(zw), save='plots/gals_xyz.png')
```


<!-- fe galaxy table -->

----------------------------------------------------------------
# View run stats: report_time calculation times
<!-- fs -->
```python
# get the df
import pandas as pd
from matplotlib import pyplot as plt
pd.set_option('display.max_columns', 30)
zrunfout='data/runtimes.dat'
zrf = pd.read_csv(zrunfout, delim_whitespace=True)

# plot mean runtimes
histcols = ['stack_boxes', 'push_box2catz', 'get_ra_dec_z', 'bin_redshifs', 'get_randoms', 'push_box2catz_Rands', 'get_ra_dec_z_Rands', 'bin_redshifs_Rands', 'calc_wtheta', 'galgal_counts', 'randrand_counts', 'galrand_counts', 'counts_to_cf']
mns = zrf[histcols].mean(axis=0)
plt.figure()
mns.plot(kind='bar')
plt.ylabel('Runtime [min]')
plt.show(block=False)

# plot histogram
plt.figure()
zrf[histcols].hist()
plt.show(block=False)

```

- does not work with new file format:
-or get df and plot at same time-

```python
import myplots as mp
zrunfout='data/zruntime.dat'
zrf = mp.getplot_zruntimes(zrunfout=zrunfout) # get zrun calc times as DF and plot
```
<!-- fe report times -->


----------------------------------------------------------------
# other:
<!-- fs -->
```python
import setup as su
if su.cosmo is None:
    su.load_cosmo() # loads global cosmo object plus H0, Om0
su.load_popmock() #
```

Use do_mock_wtheta.getmock_calcwtheta() to populate a mock, stack boxes,
transform coordinates, calculate wtheta, and write the results to a file.

<!-- fe other -->



----
# Tutorials and Notes
[Getting Started](https://halotools.readthedocs.io/en/latest/quickstart_and_tutorials/getting_started_overview.html)

Length units are comoving and assumed to be in Mpc/h throughout Halotools.

----
# The following is old:
<!-- fs -->

conda create -n halotools_env astropy numpy scipy h5py requests beautifulsoup4 cython python=3.6.8


### Run notebook
need "conda install nb_conda" to use env in notebook
source activate halotools_env
<!-- conda activate halotools_env -->
jupyter notebook
source deactivate

### To convert jupyter notebook to script
jupyter nbconvert --to python HaloTools.ipynb
<!-- first may need: conda install -c conda-forge mistune -->



## Modularized
1. Use setup_halos.py to load halo catalog (make sure it has been previously cached in the system) and HOD model.
From command line:
    screen
    source activate halotools_env
    ipython
        FROM HERE CAN RUN 'import startup' or '%run ./startup.py'
        SCRIP INCLUDES MOST OF FOLLOWING COMMANDS (BEFORE PLOTS)
    import setup_halos as sh
    halocat, HODmodel = sh.setup_halos()

2. Use calc_wtheta.py
    import calc_wtheta as cw
    import numpy as np
    make bins. e.g.:
        bins = np.logspace(np.log10(5.0), np.log10(15.0), 15)
    bcens, wtheta = cw.get_wtheta(halocat, HODmodel, bins, fout='wtheta.dat')

3. Plots
    %matplotlib qt
    import myplots as mp

    galaxy_table = HODmodel.mock.galaxy_table
    mp.plot_galaxies(galaxy_table, gal_frac=0.05, coords='radecz')

    mp.plot_wtheta(bcens, wtheta)


<!-- fe old -->
