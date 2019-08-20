# BAO

# Example Usage

## Calculate stats on 10 mocks and plot results
```python
import main as main
import plots as plots
param_dict = {  'sim_name':'multidark',
                'stats': ['wtheta', 'xi', 'wp'],
                'galplots':False}
for i in range(10):
    fstats = main.proc_mockbox(param_dict)

plots.plot_stats(fstats, save=None, show=True)
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
