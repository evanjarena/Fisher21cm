from classy import Class
from astropy.cosmology import Planck15
import numpy as np

"""8 parameter model: Some cosmological parameters are taken from 
    https://wiki.cosmos.esa.int/planckpla2015/index.php/Cosmological_Parameters
    file: https://wiki.cosmos.esa.int/planckpla2015/images/f/f7/Baseline_params_table_2015_limit68.pdf
Setion 2.7 Planck +EE + BAO, bestfit values, 
while others are obtained from astropy.cosmology Planck15
"""

class Params():

    def __init__(self,biases):
        self.biases=biases

    def param_names(self):
        """Names of cosmological parameters in LaTeX format to be used in triangle plot
        """
        names=param_names=[r"$\tau_{\mathrm{reio}}$", 
                           r"$\Omega_{\mathrm{cdm}}h^2$", 
                           r"$A_{\mathrm{s}}$",
                           r"$\theta_{\mathrm{s}}$", 
                           r"$N_{\mathrm{eff}}$",
                           r"$m_{\mathrm{ncdm}}$", 
                           r"$\Omega_{\mathrm{b}} h^2$", 
                           r"$n_{\mathrm{s}}$"]
        return names

    def default_params(self):
        """The 8 parameter cosmological model
        """
        tau_reio=0.0865            #Reionization optical depth
        omega_cdm=0.11910          #Cold dark matter density Omega_{cdm}h^2
        A_s=2.234e-9               #Amplitude of the power spectrum
        theta_s_100=1.042143       #Sound horizon at last scattering
                                   # (alternative to h)
        N_ur=Planck15.Neff-1       #Number of relativistic species in young universe
                                   # (effictive number of neutrino species)
                                   #Note this is Neff, which is now depracated
        N_ncdm=+1                  #Number of heavy neutrinos 
                                   #Needed to specify m_ncdm
        m_ncdm=0.06                #Mass of each neutrino species (in eV)
                                   #Note that this is mnu.
        omega_b=0.022319           #Baryon density Omega_{b}h^2
        n_s=0.96708                #Index of the power spectrum (spectral index)
 
        pars=np.array((tau_reio,
            omega_cdm,
            A_s,
            theta_s_100,
            N_ur,
            #N_ncdm, 
            m_ncdm,
            omega_b,
            n_s))

        pars=np.array(list(pars)+list(self.biases))
        return pars

    def perturbed_params (self,i,step):
        """Varies parameter values in order to calculate derivatives
        """
        pars=self.default_params()
        pars[i]+=step
        return pars

    def steps(self):
        """Defines step size for numerical differentiation
        """
        return self.default_params()*0.01
