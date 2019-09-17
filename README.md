# BAO

# Example Usage

## Calculate stats on N mocks and plot results (on Korriban)
```bash
printf "\e[?2004l" # turns off "bracketed paste mode" for correct copy-paste
source .bashrc
cd Documents/BAO/
mounto # needed to load Outer Rim data from Osiris
htenv
ipython
```
```python
import main as main
import plots as plots
param_dict = {  'sim_name': 'outerrim',
                'stats': ['wtheta', 'xi', 'wp'],
                'galplots': False,
                'keep_halo_frac': 0.1,
                'zbin_width': 0.05 }
for i in range(25):
    fstats = main.proc_mockbox(param_dict)

# look at stats file
p = main.load_param_dict(param_dict)
df = plots.load_statsdat(str(p['statfout']))
plots.plot_stats(str(p['statfout']), save=None, show=True, param_dict=p)
plots.plot_stats(str(p['statfout']), p, save=None, show=True, zbin='sep')

# move default file statfout
p = main.load_param_dict(param_dict)
main.file_ow(p['statfout'])
```


### Get param_dict
```python
import main as main
param_dict = {  'sim_name': 'outerrim',
                'stats': ['wtheta', 'xi', 'wp'],
                'galplots': False}
p = main.load_param_dict(param_dict)
```


## Load halo catalog
```python
import main as main
import get_sim as gs
p = main.load_param_dict({'sim_name':'multidark'})
halocat = gs.load_multidark(p)
```
