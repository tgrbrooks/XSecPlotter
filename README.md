# XSecPlotter

Cross section plotting and analysis tool for the Short-Baseline Near Detector neutrino experiment.

## Dependencies

* ROOT 6
* Input files produced by the `XSecTree` module in `sbndcode`

## Installation

The program should work on any system with standard ROOT 6 installed. Installation instruction can be found [here](https://root.cern.ch/building-root).

There are no other non standard dependencies or optional ROOT modules required.

## Running

The software is configured with the file `config.txt`, this is the only file that needs to be modified.

The code can be run either in interpreter mode with

```bash
root -l -b -q XSecPlotter.C
```

or it can be compiled with

```bash
root -l -b -q "XSecPlotter.C+"
```

It is recommended to compile the code when calculating systematics as it will significantly speed things up.

## Features

The configuration options are briefly described in the `config.txt` file.

The main features of the configuration are:
* Plot multiple cross section predictions by specifying more than one `InputFile`
* Select neutrino interaction topologies
  * Neutrino flavour
  * Charged current or neutral current interactions
  * Select by final state topology or true interaction mode
  * Define a fiducial volume
  * Select based on particle containment
* Choose which stage of reconstruction to plot with `Stage`
  * Truth level information
  * Particle reconstruction efficiencies
  * Kinematic variable smearing
  * Reconstructed selection (parametrised based on full SBND simulations) (only for numuCC)
* Choose the kinematic variable to produce differential cross sections in (supports up to 2)
* Scale to desired protons on target (POT)
* Plot rate or cross section predictions
* Histogram binning options
  * Set ranges for each parameter
  * Set number of bins or provide bin edges
  * Define a maximum statistical error per bin for automatic rebinning
* Histogram style options
  * Stack histograms by true FSI, interaction type or neutrino type
  * Show error bars on the histogram or as a percentage error band below the histogram
* Statistical analysis
  * Calculate cross section, flux, detector, external background, and constant systematic uncertainties on both rate predictions and expected cross section measurements.
  * Handle statistical uncertainty scaling with POT.
  * Calculate goodness-of-fit between models using chi2 statistical test for correlated uncertainties.
* Plotting options
  * Plot rate and cross section predictions
  * Plot 1D slices of 2D histograms
  * Plot systematic universe variations
  * Plot covariance and correlation matrices
  * Plot selection efficiency and purity graphs
  * Plot response matrices
