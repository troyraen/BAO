# Run notebook
need "conda install nb_conda" to use env in notebook
source activate halotools_env
jupyter notebook
source deactivate

# To convert jupyter notebook to script
ipython nbconvert --to python HaloTools.ipynb
first may need: conda install -c conda-forge mistune


# Tutorials
[Getting Started](https://halotools.readthedocs.io/en/latest/quickstart_and_tutorials/getting_started_overview.html)
