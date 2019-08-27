#

<!-- fs Outer Rim Setup -->
# Load Outer Rim data on Korriban
<!-- fs using Argonne GenericIO -->
```python
import pandas as pd
import genericio as gio
import plots
halo_dir = "/home/tjr63/Osiris/BAO_simdata/OuterRim/M000/L4225/HACC000/analysis/" + \
            "Halos/b0168/FOFp/"
# redshift 0.502242
halo_files = ['STEP331/m000-331.fofproperties'] + \
             ['STEP331/m000-331.fofproperties#{}'.format(n) for n in range(87,105)]
# redshift 0.539051
halo_files = ['STEP323/m000-323.fofproperties'] + \
          ['STEP323/m000-323.fofproperties#{}'.format(n) for n in range(87,105)]

f = halo_dir+halo_files[0]
gio.inspect(f)

# Load all files and ranks to a DataFrame
read_cols = [ 'fof_halo_tag', 'fof_halo_mass', \
              'fof_halo_center_x', 'fof_halo_center_y', 'fof_halo_center_z' ]
data = gio.read(f, read_cols)
df_cols = [ 'halo_id', 'halo_mass', 'x','y','z' ]
df = pd.DataFrame(data.T, columns=df_cols)
# Plot coords of a random sample (see plot below)
plots.plot_galaxies(df, title="Outer Rim z=0.539 Halos (Argonne)")
# now have ALL data loaded

# Generate halocat and populate with HOD
selct = data.T[:,2] < 500 # get halos with x<500 Mpc/h
df_cols = [ 'halo_id', 'halo_mass', 'halo_x','halo_y','halo_z' ]
df = pd.DataFrame(data.T[selct,:], columns=df_cols)
df['halo_upid'] = -1
df.drop(columns='halo_count', inplace=True)
metadata, halodata = load_outerrim_halotools_setup(p, df)
halocat = UserSuppliedHaloCatalog(**metadata, **halodata) # gives MemoryError
```
<img src="plots/outerrim/z0.539_halos_argonne.png" alt="z0.539_halos_argonne" width="600"/>
<!-- fe using Argonne GenericIO -->

<!-- Using the code below (with forked genericio) does NOT load all downloaded data. -->
<!-- fs using forked version of Generic_IO -->
```python
import sys
import os
# sys.path.insert(1, os.path.join(sys.path[0], '/home/tjr63/Documents/genericio/python/'))
# from mpi4py import MPI
import generic_io
import numpy as np
import pandas as pd
import plots # assumes in BAO dir

halo_dir = "/home/tjr63/Osiris/BAO_simdata/OuterRim/M000/L4225/HACC000/analysis/" + \
            "Halos/b0168/FOFp/"
# redshift 0.502242
halo_files = ['STEP331/m000-331.fofproperties'] + \
             ['STEP331/m000-331.fofproperties#{}'.format(n) for n in range(87,105)]
# redshift 0.539051
halo_files = ['STEP323/m000-323.fofproperties'] + \
          ['STEP323/m000-323.fofproperties#{}'.format(n) for n in range(87,105)]

# Load and look at the non-hashed file
f = os.path.join(halo_dir,halo_files[0])
gio = generic_io.Generic_IO(f, None)
metadata = gio.read_metadata(1)
colnames = gio.read_column_headers()

# Load all files
read_cols = [ 'fof_halo_count', 'fof_halo_tag', 'fof_halo_mass', \
              'fof_halo_center_x', 'fof_halo_center_y', 'fof_halo_center_z' ]
coldat = []
for hf in halo_files:
    gio = generic_io.Generic_IO(os.path.join(halo_dir,hf), None)
    coldat.append(gio.read_columns(read_cols))
df = pd.concat(coldat, axis=0)
df = df.rename(
            columns={'fof_halo_center_x':'x', 'fof_halo_center_y':'y', 'fof_halo_center_z':'z'})

# Plot coords of a random sample (see plot below)
plots.plot_galaxies(df, title="Outer Rim z=0.539 Halos")

# Inspect a random file... am I loading all the data? Seems like no.
## bash
## wc -l /home/tjr63/Osiris/BAO_simdata/OuterRim/M000/L4225/HACC000/analysis/Halos/b0168/FOFp/STEP323/m000-323.fofproperties#90
## 13404298
f = os.path.join(halo_dir,'STEP323/m000-323.fofproperties#90')
gio = generic_io.Generic_IO(f, None)
cdat = [gio.read_columns(colnames=read_cols)]
df1 = pd.concat(cdat, axis=0)
len(df1)
# 1534862
## OoM lower than number of lines in the file!
md, coords, dims = ([] for i in range(3))
for hf in halo_files:
    gio = generic_io.Generic_IO(os.path.join(halo_dir,hf), None)
    metadata = gio.read_metadata(0)
    numranks = metadata['num_ranks']
    for i in range(numranks):
        m = gio.read_metadata(i)
        coords.append(np.array(m['coords']))
        # dims.append(np.array(m['dims'])) # these are ALL [12,10,10]
        md.append(m)
cdf = pd.DataFrame(np.vstack(coords), columns=['x','y','z'])
# ddf = pd.DataFrame(np.vstack(dims), columns=['x','y','z'])
plots.plot_galaxies(cdf, title="'coords', every rank and file") # see plot below

coldat = []
for i in range(numranks):
    coldat.append(gio.read_columns(colnames=read_cols))
df = pd.concat(coldat, axis=0)
dup = df.duplicated()
sum(dup)
#######
```
<img src="plots/outerrim/z0.539_halos.png" alt="z0.539_halos" width="600"/>
<img src="plots/outerrim/coords_all_rankfile.png" alt="coords_all_rankfile" width="600"/>
<!-- fe using forked version of Generic_IO -->


<!-- fe Outer Rim Setup -->
