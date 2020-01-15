# To Do
- Understand each correlation fnc and how corrfunc calculates them
- What affects noise vs. signal strength? (zbin width, Nrands, halo downsampling)
    - try several values and make plots
- Jackknife or bootstrap to estimate errors
- Calc wp Rongpu's way (LRG clustering paper)
- consider getting rid of halos with mass less than Mmin - a few * sigma_logM
- Dark setup

# Dark setup status 9/18
i think genericio has installed correctly
conda halotools_env will not install. probably a lack of space issue
use `scl enable devtoolset-8 bash` to install with GCC 8

# Jeff 9/17
can bootstrap spacial regions (as opposed to objects)
 - only need to calculate pair counts in each region once, then combine based on bootstrap draws (with replacement)
 - try zwidth = 0.05. if still don't see signal in wtheta, can just use wp but calc this same way Rongpu does... use buffer regions to eliminate finger of god and keiser effect issues (this is what's blurring out signal in wtheta.. angle -> distance different at front of bin than back of bin.). get cross correlations (D1D2, D1R2, D2R1, R1R2) where 1 and 2 are different redshift bins.
 - try zwidth = 0.3 but with HOD_params lookup according to sim_redshift (done, no difference)

# Andrew 9/13
see if wtheta without last (last few?) zbins is less noisy (done, yes much less noisy)
estimate errors.
    - error propagation (see handwritten notes).
    - jackknife - how does covariance vary as fnc of jackknife cell size
    - see if halotools correlation functions give similar results

# Andrew 9/6
CMB and BBN tightly constrain w_b
w_b = Omega_b*h^2
uncertainty in w_b is mostly from cosmic variance, it is less than 1%
consider getting rid of halos with mass less than Mmin - a few * sigma_logM


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
scl enable devtoolset-8 bash
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


<!-- fs General HOD  -->
# Rongpu HOD params
Slack msg (1/15/19): "Hi Troy, here are the HOD parameters that I estimated for DESI LRGs"
             alpha       logM1       sigma_logM     logM0     logMmin
0.4<z<0.5    1.32828278 14.03372441  0.23939566 12.26827089 12.89376701
0.5<z<0.6    1.28471955 14.07517928  0.10847995 12.56415841 12.96064825
0.6<z<0.7    1.27655781 13.96366114  0.11947844 12.52597874 12.88283218
0.7<z<0.8    1.31175508 14.03983411  0.23792983 12.11739015 12.9599072
0.8<z<0.9    1.18976404 14.36159006  0.37759101 12.3080901  13.21897122

# Zheng07 model
To see how the following parameters are implemented, see Zheng07Cens.mean_occupation.
• param_dict[‘logMmin’] - Minimum mass required for a halo to host a central galaxy.
• param_dict[‘sigma_logM’] - Rate of transition from ⟨Ncen⟩=0⇒⟨Ncen=1⟩.
To see how the following parameters are implemented, see Zheng07Sats.mean_occupation.
• param_dict[‘alpha’] - Power law slope of the relation between halo mass and ⟨Nsat⟩.
• param_dict[‘logM0’] - Low-mass cutoff in ⟨Nsat⟩.
• param_dict[‘logM1’] - Characteristic halo mass where ⟨Nsat⟩ begins to assume a power law form.
<!-- fe General HOD -->

<!-- fs Conda Environment Setup -->
# Setup Conda Environment
```bash
# scl enable devtoolset-8 bash
conda create -n halotools_env astropy numpy scipy h5py requests beautifulsoup4 cython python=3.7
conda activate halotools_env
pip install halotools
# install Generic_IO (see above)
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
cd ${OGdir}
```
<!-- fe Conda Environment Setup -->


# Dark
Scientific Linux 6.10
ssh <username>@dark.phyast.pitt.edu
