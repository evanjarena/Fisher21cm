# Fisher21cm

| Fisher21cm |-----| Fisher forecast for a general 21cm experiment |
| Authors    |-----| Evan J. Arena and Anze Slosar                  |

* (c) Evan J. Arena (Drexel University Department of Physics), 2020.
* For questions, please email `evan.james.arena@drexel.edu.`

## Required Packages

* numpy
* scipy
* matplotlib
* h5py
* pickle
* SNR data cube from Cosmic Visions (optional, see below)
* CLASS and/or CAMB Boltzmann codes:
   * [CLASS code and its python wrapper classy](http://class-code.net/)
   * [CAMB code and its python wrapper](http://camb.info)

## Modules

Modules should be run from the parent directory, e.g.:
  `python py/21cm.py Class`
or
  `python py/21cm.py Camb`
in order to run the Fisher forecast with either the CLASS or CAMB Boltzmann codes.

* `21cm.py`: Details for a generalized 21cm experiment.  The default is a BMX-type experiment, but the user can replace this module with any general experiemnt. 

* `ParameterVec.py`: Holds the cosmological parameters.

* `FisherMat.py`: Fisher Matrix class.

* `TracerPk.py`: Matrix for a tracer measuring the matter power spectrum.

* `FisherAnalysis.py` - Computes the covariance matrix and creates triangle plots of 2D marginalized error ellipses and 1D marginalized Gaussians.  Can combine Fisher matrices to compare experiments.

* `S4Fisher.py`: Gets CMB-S4 Fisher matrix for comparison to 21cm experiments.  

* `ClassWrap.py`: Wrapper for the Boltzmann code CLASS.

* `CAMBWrap.py`: Wrapper for the Boltzmann code CAMB.

By default, the signal-to-noise ratio (SNR) for the Fisher Matrix is computed in a simulation from Cosmic Visions:
* [SNR Data Cube from Cosmic Visions](http://www.phas.ubc.ca/~richard/sn_lowz_expA_50K.h5)
* [Cosmic Visions Fisher repo](https://github.com/radiohep/CVFisher)



