# PySimSpin

PySimSpin v2.0.0 - A package for producing mock observations of particle simulations

The purpose of the PySimSpin python package is to take a particle simulation of a galaxy and produce a spectral data cube in the style of a specified Integral Field Spectroscopy (IFS) instrument. A mock spectral data cube can be produced using the functions in this package. This is a simple process comprised of three steps:

1. Read in your particle data and produce the relevant spectra using the `make_simspin_file` function.
1. Setup the observation by defining your `telescope` and `observing_strategy`.
1. Build your data cube using the `build_datacube`.

From this cube, "observables" can be measured using observational pipelines. Specifically, the observed line-of-sight velocities and dispersions. This package, once installed, is fully documented and tested.

Another implementation of this code (SimSpin v1.1.3) written in Julia is also available at [SimSpin.jl](https://github.com/kateharborne/SimSpin.jl) developed by [Gerry Gralton](https://github.com/gerrygralton).

A further R-package implementation of SimSpin v2.0.0 is also in development at [SimSpin](https://github.com/kateharborne/SimSpin).
