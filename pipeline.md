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


# Modularized
1. Use setup_halos.py to load halo catalog (make sure it has been previously cached in the system) and HOD model.
From command line:
    source activate halotools_env
    ipython
    import setup_halos as sh
    halocat, HODmodel = sh.setup_halos()

2. Use get_stats.py
    import get_stats as gs
    make bins. e.g.:
        bins = np.logspace(5., 15., num=10)
    gs.get_stats(halocat, HODmodel, bins)
