# Run notebook
need "conda install nb_conda" to use env in notebook
source activate halotools_env
jupyter notebook
source deactivate

# To convert jupyter notebook to script
jupyter nbconvert --to python HaloTools.ipynb
<!-- first may need: conda install -c conda-forge mistune -->


# Tutorials
[Getting Started](https://halotools.readthedocs.io/en/latest/quickstart_and_tutorials/getting_started_overview.html)

----
# Modularized
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
