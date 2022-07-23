# drft

drft is a package to modify
[RAMSES](https://bitbucket.org/rteyssie/ramses) (Teyssier 2002)
initial conditions (ICs) made by
[MUSIC](https://bitbucket.org/ohahn/music/src) (Hahn & Abel 2011) to
include the effects of relative baryon-dark matter velocities
(Tseliakhovich & Hirata 2010) from recombination to the starting time
of your simulation. It contains routines for solving the evolution
equations (`py_vbc`), which is essentially a Python reproduction of
`vbc_transfer` in CICsASS (O'Leary & McQuinn 2012). drft was
itself based on routines in
[seren3](https://github.com/sully90/seren3) written by David Sullivan.


### How to use

The basic steps for running drft are:

1. Compute transfer functions that have separate amplitudes for baryon
and dark matter velocities (e.g. using
[CAMB](https://github.com/cmbant/CAMB/), Lewis 2000)

2. Use these transfer functions to generate cubic ICs using [MUSIC](https://bitbucket.org/ohahn/music/src)
(this can be enforced by setting `force_equal_extent=yes` in the
[MUSIC](https://bitbucket.org/ohahn/music/src) configuration file)

3. Calculate the v_bc field (this can be done using `yt_ic_vbc.py`,
which uses [yt](https://github.com/yt-project/yt), Turk et al. 2011)

4. Run `bias_ics.py` (see `work_ics.sh.example` for an example script)

Cosmological parameters for `py_vbc` can be set using a `.ini` file,
see `py_vbc/planck2018_params.ini` for an example. If you do want to
change parameters, you will need to:

1. Change `config_fname` in `py_vbc/constants.py` to the new
configuration filename

2. Regenerate transfer functions at z=1003, 1000, 997 and 200 (for the
setup described in the paper) and place them in `py_vbc/tfs/`

3. Regenerate a
[RECFAST](https://www.astro.ubc.ca/people/scott/recfast.html)
(Seager et al. 1999) output and place it in `py_vbc/recfast/`


### Acknowledging

If you use results produced by this package in a scientific
publiocation, please cite the methodology paper.


### References

Hahn O., Abel T., 2011, Monthly Notices of the Royal Astronomical
Society, 415, 2101

Lewis A., Challinor A., Lasenby A., 2000, The Astrophysical Journal,
538, 473

O’Leary R. M., McQuinn M., 2012, The Astrophysical Journal, 760, 4

Seager S., Sasselov D. D., Scott D., 1999, The Astrophysical Journal
Letters, 523, L1

Teyssier R., 2002, Astronomy and Astrophysics, 385, 337

Tseliakhovich D., Hirata C., 2010, Physical Review D, 82, 083520

Turk M. J., Smith B. D., Oishi J. S., Skory S., Skillman S. W., Abel T., Norman
M. L., 2011, The Astrophysical Journal Supplement Series, 192, 9