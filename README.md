# Stellar population synthesis and stellar mass fitting with BC03 models

This module uses the stellar population synthesis code GALAXEV (Bruzual \& Charlot 2003, BC03) to create composite stellar population (CSP) models and predict broad band magnitudes on a grid of stellar population parameters (stellar mass, star formation history, metallicity, dust attenuation). 
It also contains example scripts to fit for the stellar mass of galaxies observed with VST, using BC03 CSP models.

## Requirements

- Python
- emcee (for the fitting part)
- bc03

To run the full set of examples, you'll need to define the environmental variable PYGALAXEVDIR, pointing to this directory.

## BC03

bc03 consists of:
- synthetic spectra corresponding to simple stellar population (i.e. single burst, SSP), obtained with either empirical (MILES, Stelib) or theoretical (BaSeL) stellar spectral libraries, and stellar evolution models.
- code that combines SSPs to create composite stellar population models with custom star formation history.

You can download bc03 from [here](http://www.bruzual.org/bc03/). Several versions of the code are available, but this module has been written to work with the 2016 version. You will need to download and unpack the 'src.tgz' and at least one of the SSP libraries (for instance, 'BC03_basel_kroupa.tgz'). Untar everything and make sure to unzip the `*.ised` files in the `Chabrier_IMF/` (or Kroupa or Salpeter) directory.

### Setting up BC03

1. Define a `$bc03` environmental variable pointing to the `src/` directory of bc03.
2. Edit the file `$bc03/.bc_bash`, changing the line containing `export bc03` to match your installation of bc03 (I suppose deleting that line would also work, if you've done point 1. above). If using `.bc_cshrc` instead, you should probably comment out line 8, `setenv STILTS`.
3. `cd $bc03`; `source ./.bc_bash`; `make all`.

In order to be able to run bc03 code properly from any directory, you'll need to source the `.bc_bash` file. You might want to add `source $bc03/.bc_bash` to your `.bashrc` file.

### Example 1: create a CSP model with exponentially decaying star formation history

The file `examples/make_one_csp.py` included in this module reads a set of SSP models from the bc03 package (in this case the models obtained with the BaSeL spectral library and a Chabrier IMF) and creates a composite stellar population model with an exponentially decaying star formation history, a given metallicity (out of the values sampled by the Padova 1994 models), a fixed exponential decay timescale 'tau' and a fixed dust optical depth tau_V.
To run it,
- add this directory to your PYTHONPATH environment variable
- copy the file to a directory of your choice
- change the path assigned to the variable `ssp_dir` in line 8, to match that of your installation of bc03 and your desired SSP models.
- change the path assigned to the variable `work_dir` in line 9, to wherever you wish to run this code.
- run `python make_one_csp.py`

The code should take a few seconds to run. The end result should be a pair of files named 'bc03_chab_...eps=0.000.mass' and 'bc03_chab_...eps=0.000.ised'. The `.ised` files contain the synthetic spectral energy distribution of the CSP over a broad range of ages, but are in binary format. To convert them to a more readable format, you'll need to run the script `galaxevpl` included in bc03 (but the code presented in the next example will do that for you). The `.mass` files are tables listing the evolution of the stellar mass as a function of age for each CSP. They also list luminosity and M/L in a few standard bands.

### Example 2: store model SED in .hdf5 files

The code `examples/get_one_sed_from_csp.py` reads in the files produced by running the `make_one_csp.py` script, renormalizes the model SEDs to a total stellar mass (living stars + remnants) of 1 solar mass, and stores it into .hdf5 files. As in the previous example, copy this file to your preferred directory and then update the values of `ssp_dir` and `work_dir`. This code should run very quickly. The end product should be a new files named `bc03_Z=....age=13.000.hdf5`. This is a synthetic spectrum. The file size depends on the resolution of the initial SSP models (see note at the bottom on the choice of SSP models).

### Example 3: calculate broad band magnitudes of a galaxy given a model CSP SED.

Now that we have a synthetic spectrum for our CSP, we can use it to predict model magnitudes of a given galaxy. The script `examples/get_mags_from_sed.py` calculates the magnitude in VST u, g, r, i bands for a galaxy at a given redshift and stellar mass and prints them to the standard output.
**The code assumes flat LCDM with OmegaM=0.3 and H0=70 by default**.

### Example 4: creating a grid of CSP SEDs.

It can be convenient to produce synthetic spectra over a grid of stellar population parameters in a single step, for later use. This is done by the script `examples/make_sed_grid.py`. This code takes hours to run. The main output will be a .hdf5 file containing all spectra at each point of the grid. The script also outputs a large list of files of the kind of those produced by Example 1 (one for each grid point in metallicity, tau, and tau_V). These may be deleted.

### Example 5: calculate broad band magnitudes on a grid of stellar population parameters.

The script `examples/make_mags_grid_fixedz.py` calculates broad band magnitudes (for a 1 Solar mass galaxy) in a set of filters at each point of the stellar population parameter grid created in example 4, at a fixed redshift. It should take a few minutes to run.

### Example 6: fit CSP models to the observed broad band magnitudes to infer the stellar mass and the other stellar population parameters
The script `examples/fit_bbmags.py` runs an MCMC to fit for the SPS parameters of a galaxy, including the stellar mass, using the grid of model magnitudes created with the `examples/make_mags_grid_fixedz.py` script. The output is an .hdf5 file with the chain, as produced by emcee.

### Notes on the choice of SSP models

The 2016 version of bc03 comes with three sets of SSP, obtained with different stellar libraries. The 'BaSeL' libraries are low-resolution (20\AA) theoretical spectra, while the 'Stelib' and 'MILES' are based on observed high-resolution (3\AA) spectra. Since this code is meant to be used to fit broad-band photometry data, using high-resolution spectra is overkill and results in much larger file sizes and computational time. So, if you wish to use either the 'Stelib' or 'MILES' SSPs, I recommend first downgrading the resolution by rebinning the spectra by running the script `downgrade_resolution` (see section 3.6 in the [bc03 documentation](http://www.bruzual.org/bc03/doc/bc03.pdf) for details).

### Notes on photometric uncertainties

KiDS can measure mangitudes with very high precision (formal uncertainties smaller than 0.01). However, the CSP models are not as accurate. So, if you try to use this code to fit KiDS magnitudes with uncertainties taken directly from the catalog, what typically ends up happening is that no model describing the data can be found. The MCMC chain will oscillate around a point with very small fluctuations, as is usually the case when combining highly inconsistent measurements (i.e. the product of two Gaussian distribution with means that are several sigmas apart from each other is a very narrow Gaussian).
I therefore recommend adding a systematic uncertainty to the observed magnitudes, to account for this model mismatch. In my experience, 0.05 magnitudes should do the job.
