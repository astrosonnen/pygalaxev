# Stellar population synthesis and stellar mass fitting with BC03 models

This module uses the stellar population synthesis code GALAXEV (Bruzual \& Charlot 2003, BC03) to create composite stellar population (CSP) models and predict broad band magnitudes on a grid of stellar population parameters (stellar mass, star formation history, metallicity, dust attenuation). 
It also contains example scripts to fit for the stellar mass of galaxies observed with KiDS, using BC03 CSP models.

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

### Example 1: create CSP models with exponentially decaying star formation history

The file `examples/make_csp.py` included in this module reads a set of SSP models from the bc03 package (in this case the models obtained with the BaSeL spectral library and a Kroupa IMF) and creates composite stellar population models corresponding to an exponentially decaying star formation history, over a grid of values of metallicity, star formation decay time and dust attenuation.
To run it,
- copy the file to a directory of your choice
- change the path assigned to the variable `ssp_dir` in line 8, to match that of your installation of bc03 and your desired SSP models.
- change the path assigned to the variable `work_dir` in line 9, to wherever you wish to run this code (keep in mind that around 10k files will be created).
- add `stellarpop` to your `$PYTHONPATH` environment variable
- run `python make_csp.py`

This code should take a few hours to run. The end result should be a large (i.e. thousands) set of files named 'bc03_chab_...eps=0.000.mass' and 'bc03_chab_...eps=0.000.ised'. The `.ised` files contain the synthetic spectral energy distribution of each CSP, but are in binary format. To convert them to a more readable format, you'll need to run the script `galaxevpl` included in bc03 (but the code presented in the next example will do that for you). The `.mass` files are tables listing the evolution of the stellar mass as a function of age for each CSP. They also list luminosity and M/L in a few standard bands.

### Example 2: store model SEDs in python-readable files

The code `examples/make_sed_from_csp.py` reads in the files produced by running the `make_csp.py` scripts, renormalizes the model SEDs to a total stellar mass (living stars + remnants) of 1 solar mass, and stores them into 'pickled' files. As for the previous example, copy this file to your preferred directory and then update the values of `ssp_dir` and `work_dir`. This code should take two/three hours to run. The end product should be a few files named `chabrier_Z=....dat`. These are essentially synthetic spectra, so the file size depends on the resolution of the initial SSP models (see note at the bottom on the choice of SSP models).

### Example 3: calculate broad band magnitudes on a grid of stellar population parameters for a galaxy at fixed redshift.

Now that we have synthetic spectra for our CSPs, we can use them to predict model magnitudes for a given galaxy. The script `examples/make_models_fixedz.py` creates an object that allows one to evaluate the magnitude in KiDS u, g, r, i and VISTA k bands for a galaxy at fixed redshift, by interpolating over a grid. This is done by running the script `examples/get_mags_from_model.py`, which calculates the model magnitudes of a galaxy at z=0.609 with given stellar mass, age, star formation decay time, metallicity, dust attenuation. 
**The code assumes flat LCDM with OmegaM=0.3 and H0=70 by default**.

### Example 4: fit for the stellar population parameters of a given galaxy

The code `examples/fit_mags_fixedz_model.py` uses the model created in the previous example to fit for the stellar mass and other parameters of a galaxy observed with HSC. As usual, please update the path in line 14 before running the script. This script requires pymc to run. It should take about a minute to run. The output MCMC chain is stored in a 'pickled' file. The chain has been 'thinned' (only 1 in 10 samples are kept, the others are discarded) to keep file size under control.

### Example 5: create a grid of models over a redshift range, to fit a large number of objects

Examples 3 and 4 are fine if you're only interested in fitting a few galaxies. If you're dealing with hundreds of objects or more, it is more efficient to add the 'redshift' dimension to the grid of models used to compute magnitudes.
The script `examples/make_models_zgrid.py` creates such a model. This should take a couple of hours to complete.
When it's done, you can run the script `examples/get_mags_from_zgrid_model.py` to evaluate model magnitudes for a galaxy at any redshift between z=0.01 and z=1.0.

### Example 6: fit for the stellar population parameters of a given galaxy, using the model grid created in step 5
The script `examples/fit_mags_zgrid_model.py` runs an MCMC to fit for the SPS parameters of a galaxy, using the grid of model magnitudes created with the `examples/make_models_zgrid.py` script.

### Notes on the choice of SSP models

The 2016 version of bc03 comes with three sets of SSP, obtained with different stellar libraries. The 'BaSeL' libraries are low-resolution (20\AA) theoretical spectra, while the 'Stelib' and 'MILES' are based on observed high-resolution (3\AA) spectra. Since this code is meant to be used to fit broad-band photometry data, using high-resolution spectra is overkill and results in much larger file sizes and computational time. So, if you wish to use either the 'Stelib' or 'MILES' SSPs, I recommend first downgrading the resolution by rebinning the spectra by running the script `downgrade_resolution` (see section 3.6 in the [bc03 documentation](http://www.bruzual.org/bc03/doc/bc03.pdf) for details).

### Notes on photometric uncertainties

KiDS can measure mangitudes with very high precision (formal uncertainties smaller than 0.01). However, the CSP models are not as accurate. So, if you try to use this code to fit KiDS magnitudes with uncertainties taken directly from the catalog, what typically ends up happening is that no model describing the data can be found. The MCMC chain will oscillate around a point with very small fluctuations, as is usually the case when combining highly inconsistent measurements (i.e. the product of two Gaussian distribution with means that are several sigmas apart from each other is a very narrow Gaussian).
I therefore recommend adding a systematic uncertainty to the observed magnitudes, to account for this model mismatch. In my experience, 0.05 magnitudes should do the job.
