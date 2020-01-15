# BAO

# Project Overview

The expansion history of the universe can be measured using the Baryon Acoustic Oscillation (BAO) signal. Precise measurements will probably require spectroscopic data from large volumes of the sky. However, spectroscopy is time-intensive and therefore large-scale surveys (like the upcoming LSST) rely on photometry which contains less information. This project (currently in the beginning stages) explores the use of photometric data and machine learning techniques to predict where on the sky spectroscopy is likely to measure the BAO with high signal-to-noise. In turn, this can be used to inform spectroscopic survey strategies so that time and resources are not wasted on regions unlikely to produce useful data.


------------------------------
# Example Usage

## Calculate stats on N mocks and plot results (on Korriban)

```bash
printf "\e[?2004l" # turns off "bracketed paste mode" for correct copy-paste
cd Documents/BAO/
mounto # needed to load Outer Rim data from Osiris
screen
conda activate halotools_env
ipython
```


```python
import main as main
import plots as plots
from pathlib import Path
param_dict = {  'sim_name': 'outerrim',
                'stats': ['wtheta', 'xi', 'wp'],
                'keep_halo_frac': 0.1,
                'zbin_width': 0.05,
                'statfout': Path().resolve() / 'data/stats_zw005.dat'
                }
N = 100
for i in range(N):
    fstats = main.proc_mockbox(param_dict)

# look at stats file, plot results
p = main.load_param_dict(param_dict)
df, dfmeta = plots.load_statsdat(str(p['statfout']))
plots.plot_stats(str(p['statfout']), cosmo_in=p['cosmo'], save=None, show=True, zbin='avg')

for i in range(10):
    plots.plot_stats(str(p['statfout']), cosmo_in=p['cosmo'], zbin='avg', keep_zbin=[5*i+j for j in range(5)])
plots.plot_stats(str(p['statfout']), cosmo_in=p['cosmo'], zbin='avg', keep_zbin=[5*10+j for j in range(5)])

plots.plot_stats(str(p['statfout']), cosmo_in=p['cosmo'], zbin='avg', keep_zbin=[i for i in range(10,40)
```


## Some useful individual tasks

```python
# move default file statfout
p = main.load_param_dict(param_dict)
main.file_ow(p['statfout'], suffix='_zw01')

# Get param_dict with unset params filled with defaults
import main as main
param_dict = {  'sim_name': 'outerrim',
                'stats': ['wtheta', 'xi', 'wp'],
                'galplots': False}
p = main.load_param_dict(param_dict)


# Load halo catalog
import main as main
import get_sim as gs
p = main.load_param_dict({'sim_name':'multidark'})
halocat = gs.load_multidark(p)
```
