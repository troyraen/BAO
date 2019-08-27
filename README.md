# BAO

# Example Usage

## Calculate stats on 10 mocks and plot results (on Korriban)
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
                'sim_redshift': 0.539051,
                'sim_Lbox': 3000.0,
                'sim_particle_mass': 1.85e9,
                'stats': ['wtheta', 'xi', 'wp'],
                'galplots': False}
for i in range(10):
    fstats = main.proc_mockbox(param_dict)

plots.plot_stats(str(fstats), save=None, show=True)

# move default file statfout
p = main.load_param_dict(param_dict)
main.file_ow(p['statfout'])
```


### Get param_dict
```python
import main as main
p = main.load_param_dict({})
```


## Load halo catalog
```python
import main as main
import get_sim as gs
p = main.load_param_dict({'sim_name':'multidark'})
halocat = gs.load_multidark(p)
```
