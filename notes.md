# SETUP
Korriban:
```bash
source .bashrc
cd Documents/BAO
htenv
# source activate halotools_env
```


# Outer Rim Simulation Info
Downloading data with Globus Connect Personal
Korriban setup key c41f2d31-aab6-4e02-abfe-a075554cd226 # could not complete setup here
Roy setup key f6be8bf5-c3e7-47e4-85c8-5387b1980459

h = 0.71
Volume = (4.225 Gpc)^3 = (3000 h^-1 Mpc)^3
Num Particles = 10240^3
Mass resolution = 2.6e9 M_sun = 1.85e9 h^-1 M_sun
FoF linking length, b = 0.168


# General Notes:
Expect BAO at ~4.5 degrees for z~0.5 (see plots/zthetaBAO.png)
Halotools assumes all lengths are in Mpc/h
max redshift error in SDSS DR10 Photoz table is 0.365106,
    see http://skyserver.sdss.org/CasJobs/MyDB.aspx MyTable_1 and
    http://skyserver.sdss.org/dr6/en/help/docs/algorithm.asp?key=photoz
    Note from Jeff: SDSS has unusually large z errors.
