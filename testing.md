# Test NR/NG = 0.1, 0.01
<!-- fs -->
'''python -u main.py >> mainNRNGsmall.out'''
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



# Test Corrfunc.theory xi, wp
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
