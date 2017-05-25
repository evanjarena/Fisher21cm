# Fisher21cm

This is the repo for the Baryon Mapping eXperiment (BMX) Fisher forecast for a general 21cm experiment.

## Modules

`parameters.py` - The cosmological parameters in our eight parameter model.

`power_spectrum.py` - Computes the matter power spectrum 

\P(k_{\parallel},k_{\perp})=b^2(1+\beta\mu^2)^2P_L(k_T),

where \beta=f/b, f is the growth rate, b is the bias, and \mu=\cos\theta=\cos(k_{\parallel}/k_T).

`fisher_calc.py` - Computes the Fisher matrix

F_{ij}=\sum_{k_{\parallel},k_{\perp},z}\frac{\partial P(k_{\parallel},k_{\perp})}{\partial \theta_i}\frac{\partial P(k_{\parallel},k_{\perp})}{\partial \theta_j}\frac{1}{(\Delta P)^2}\\
 = \sum_{k_{\parallel},k_{\perp},z}\frac{\partial \ln P(k_{\parallel},k_{\perp})}{\partial \theta_i}\frac{\partial \ln P(k_{\parallel},k_{\perp})}{\partial \theta_j}(\text{SNR})^2.\\

Here, the signal-to-noise ratio (SNR) is computed in a simulation from Cosmic Visions:
- [SNR Data Cube from Cosmic Visions](http://www.phas.ubc.ca/~richard/sn_lowz_expA_50K.h5)
- [Cosmic Visions Fisher repo](https://github.com/radiohep/CVFisher)

`fisher_analysis.py` - Computes the covariance matrix and creates triangle plots of 2D marginalized error ellipses and 1D marginalized Gaussians.  Can combine Fisher matrices to compare experiments.


`fisherS4_pickler.py` Takes the CMBS4 Fisher matrix given in 'fisherS4.txt' and creates the pickle file 'fisherS4_matrix.p'.  This is used in fisher_analysis.py.

## Requirements
- numpy
- scipy
- matplotlib
- SNR data cube mentioned above and h5py to read it into module
- pickle
- [CLASS code and its python wrapper classy](http://class-code.net/)
