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
[CAMB](https://github.com/cmbant/CAMB/), Lewis et al. 2000)

2. Use these transfer functions to generate cubic ICs using [MUSIC](https://bitbucket.org/ohahn/music/src)
(this can be enforced by setting `force_equal_extent=yes` in the
[MUSIC](https://bitbucket.org/ohahn/music/src) configuration file)

3. Calculate the v_bc field (this can be done using `yt_ic_vbc.py`,
which uses [yt](https://github.com/yt-project/yt), Turk et al. 2011)

4. Run `bias_ics.py` (see `work_ics.sh.example` for an example script)

`py_vbc` is configured at runtime with a YAML file. A flat cosmology
is assumed. See `py_vbc/planck2018_params.yaml` for the complete
Planck example (the mapping below shows the entries needed for the
standard z=1000 run):

```yaml
cosmology:
  h: 0.673
  omega_m: 0.314
  omega_b: 0.049
  sigma_8: 0.812
  ns: 0.965
  costh: 1.0

filenames:
  transfer_functions:
    0: tfs/planck2018_transfer_out_z000.dat
    997: tfs/planck2018_transfer_out_z997.dat
    1000: tfs/planck2018_transfer_out_z1000.dat
    1003: tfs/planck2018_transfer_out_z1003.dat
```

All paths in `filenames` are resolved relative to the YAML file. Transfer
functions are explicit redshift-to-file entries, so their filenames do not
need to follow a shared prefix. For example, a run with `zstart=1000` and
`dz=3` needs entries at z=997, 1000, and 1003; power-spectrum normalization
also uses the configured z=0 entry. `rf_base` is optional; when omitted,
py_vbc uses the bundled `py_vbc/recfast/planck2018_recfast.dat` file and
prints one runtime notice. A custom `rf_base` may point to any
[RECFAST](https://www.astro.ubc.ca/people/scott/recfast.html)
(Seager et al. 1999) output. The radiation density is derived internally as
`omega_r = 4.15e-5 / (h ** 2.0)` (Dodelson 2002, Eq. 2.86), so it is not a
YAML input.

Pass the selected configuration to every run:

```python
import py_vbc

k, (p_c, p_b, p_vc, p_vb) = py_vbc.run_pyvbc(
    vbc=30.0,
    zstart=1000.0,
    zend=200.0,
    dz=3.0,
    config_file="path/to/params.yaml",
)
```

For repeated calculations, you canpreload and reuse the immutable
configuration:

```python
config = py_vbc.load_config("path/to/params.yaml")
k, power = py_vbc.run_pyvbc(..., config_file=config)
```

The main `drft` routines (`bias_ics.py`, `gen_ics.py`, and
`utils.compute_bias()`) accept the same YAML configuration; see
`work_ics.sh.example` for the command-line form used by `bias_ics.py`.


### Acknowledging

If you use results produced by this package in a scientific
publication, please cite the methodology paper Conaboy et al. (2023) [ADS](https://ui.adsabs.harvard.edu/abs/2023MNRAS.525.5479C/abstract).


### References

Conaboy, L., Iliev, I. T., Fialkov, A., Dixon, K. L., 
Sullivan, D., 2023, Monthly Notices of the Royal Astronomical
Society, 525, 5479

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
