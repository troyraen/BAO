# Initiate a Korriban session
```bash
printf "\e[?2004l" # turns off "bracketed paste mode" for correct copy-paste
source .bashrc
mounto
cd Documents/BAO
htenv
ipython
# source activate halotools_env
```

<!-- fs Outer Rim -->
# Outer Rim Simulation Info
h = 0.71
Volume = (4.225 Gpc)^3 = (3000 h^-1 Mpc)^3
Num Particles = 10240^3
Mass resolution = 2.6e9 M_sun = 1.85e9 h^-1 M_sun
FoF linking length, b = 0.168
Lightcone redshift range = [3.036145, 0.0]
Cosmology params: { Omega_CDM: 0.2200,
                    w_b: 0.02258,
                    w_nu: 0.0,
                    h: 0.7100,
                    omega_8: 0.8000,
                    n_s: 0.9630,
                    w_0: -1.0000,
                    w_a: 0.0000
                    }
<!-- fs Columns (i=integer,f=floating point, number bits size) -->
FoF columns:
    - [i 32] fof_halo_count
    - [i 64] fof_halo_tag
    - [f 32] fof_halo_mass
    - [f 32] fof_halo_center_x
    - [f 32] fof_halo_center_y
    - [f 32] fof_halo_center_z
    - [f 32] fof_halo_mean_x
    - [f 32] fof_halo_mean_y
    - [f 32] fof_halo_mean_z
    - [f 32] fof_halo_mean_vx
    - [f 32] fof_halo_mean_vy
    - [f 32] fof_halo_mean_vz
    - [f 32] fof_halo_vel_disp
Particle columns:
    - [f 32] x
    - [f 32] y
    - [f 32] z
    - [f 32] vx
    - [f 32] vy
    - [f 32] vz
    - [f 32] phi
    - [i 64] id
    - [i 16] mask
Lightcone Halos columns:
    - [f 32] x
    - [f 32] y
    - [f 32] z
    - [f 32] vx
    - [f 32] vy
    - [f 32] vz
    - [f 32] a
    - [f 32] mass
<!-- fe Columns -->

## Download data
- Install Globus Connect Personal
    * To install on Osiris, used (install GCP)[https://docs.globus.org/how-to/globus-connect-personal-linux/] instructions under 'How To Install Globus Connect Personal for Linux Using the Command Line'
    * Setup keys:
        - Korriban setup key c41f2d31-aab6-4e02-abfe-a075554cd226 # could not complete setup here
        - Roy setup key f6be8bf5-c3e7-47e4-85c8-5387b1980459
        - Osiris setup key 55a73d27-c963-4076-aedf-59a3c379a55a

- Start Globus on Osiris
    * cd /home/tjr63/Globus/globusconnectpersonal-2.3.8
    * ./globusconnectpersonal -start &
- Initiate download
    * (Argonne portal)[https://cosmology.alcf.anl.gov/outerrim]
    * Choose data to download
    * At Globus web portal:
        * 'Collection' -> 'Your Collections' -> Osiris
        * for destination dir chose 'BAO_simdata'
- Check status on Osiris
    * ./globusconnectpersonal -status
- [View web activity](https://app.globus.org/activity)


## Setup Generic_IO (to read data) in htenv on Korriban
[Argonne Generic_IO Project](https://trac.alcf.anl.gov/projects/genericio)
```bash
wget http://www.mcs.anl.gov/~turam/genericio/genericio-20190417.tar.gz
tar xvzf genericio-20190417.tar.gz
cd genericio
make frontend-progs
# add to PYTHONPATH in Conda env (see below)
```

[Forked repo](https://github.com/Christopher-Bradshaw/genericio)
```bash
htenv
git clone https://github.com/Christopher-Bradshaw/genericio.git genericio_fork # fork. instructions in readme
cd genericio_fork
make py_deps
make py_build
make py_test
# add to PYTHONPATH in Conda env (see below, but use dir genericio_fork)
```
<!-- fe Outer Rim -->


<!-- fs General BAO -->
# General BAO Notes:
Expect BAO at ~4.5 degrees for z~0.5 (see plots/zthetaBAO.png)
Halotools assumes all lengths are in Mpc/h
max redshift error in SDSS DR10 Photoz table is 0.365106,
    see http://skyserver.sdss.org/CasJobs/MyDB.aspx MyTable_1 and
    http://skyserver.sdss.org/dr6/en/help/docs/algorithm.asp?key=photoz
    Note from Jeff: SDSS has unusually large z errors.
<!-- fe General BAO -->


<!-- fs Conda Environment Setup -->
# Setup Conda Environment
```bash
conda create -n halotools_env python=3.7
conda activate halotools_env
OGdir=${PWD}
cd $CONDA_PREFIX
mkdir -p ./etc/conda/activate.d
mkdir -p ./etc/conda/deactivate.d
touch ./etc/conda/activate.d/env_vars.sh
touch ./etc/conda/deactivate.d/env_vars.sh
echo 'export OLD_PYTHONPATH="${PYTHONPATH}"' >> ./etc/conda/activate.d/env_vars.sh
echo 'export PYTHONPATH="${PYTHONPATH}:/home/tjr63/Documents/genericio/python"' >> ./etc/conda/activate.d/env_vars.sh
echo 'export PYTHONPATH=${OLD_PYTHONPATH}' >> ./etc/conda/deactivate.d/env_vars.sh
echo 'unset OLD_PYTHONPATH' >> ./etc/conda/deactivate.d/env_vars.sh
```
<!-- fe Conda Environment Setup -->
