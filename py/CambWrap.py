#
# Wrapper for CAMB.
#

from __future__ import division, print_function
import camb
from camb import model, initialpower
import copy
from scipy.interpolate import interp1d
from ParameterVec import DefaultParamList, ParamList, Parameter
import sys
import numpy as np

class PkDiffer_Camb:

    def __init__ (self,pl, zvals, kvals, kperp, kpar, Nkmu2_row, Nkmu2_col):
        """ returns a list of Pks, each list containins 3D PS"""
        #self.zstr=",".join(map(str,zvals+[zvals[-1]+2]))
        self.kvals=kvals
        self.kperp=kperp
        self.kpar=kpar
        self.zvals=zvals
        self.Nkmu2_row=Nkmu2_row
        self.Nkmu2_col=Nkmu2_col
        self.plist=copy.deepcopy(pl)
        self.pars = camb.CAMBparams()
        self.ComputeCosmo(pl)
        bg=camb.get_background(self.pars)
        nz = 4620 #number of steps to use for the radial/redshift integration
        kmax=10  #kmax to use
        #For Limber result, want integration over \chi (comoving radial distance), from 0 to chi_*.
        #so get background results to find chistar, set up arrage in chi, and calculate corresponding redshifts
        chistar = bg.conformal_time(0)- model.tau_maxvis.value
        chis = np.linspace(0,chistar,nz)
        zs=bg.redshift_at_comoving_radial_distance(chis)
        #zs=zs[::-1]
        Da=interp1d(zs,bg.comoving_radial_distance(zs)) 
        Hi=interp1d(zs,1./(bg.h_of_z(zs)))  
        self.Da_fid=Da
        self.Hi_fid=Hi
        self.cube0=self.getCube(pl,'store_fid')

    def getDerivative(self, pa, frac):
        de=[]
        # if pa.name=='theta':
        #     ufrac=frac/50.
        # else:
        #     ufrac=frac
        ufrac=frac
        npl=copy.deepcopy(self.plist)
        if (pa.value==0):
            if 'Mkmu2_' in pa.name:
                step=0.01
            elif 'Akmu2_' in pa.name:
                step=100
            else:
                print('WTF?')
                exit()
        else:
            step=pa.value*ufrac
        for fa in [+1,-1]:
            nval=pa.value+fa*step
            npl.setValue(pa.name,nval)
            if "b_" in pa.name:
                mode="use_fid"
            else:
                mode="normal"
                self.ComputeCosmo(npl)
            de.append(self.getCube(npl,mode))
        toret=[(dp-dm)/(2*step) for dp,dm in zip(de[0],de[1])]
        return toret
        
    def getCube(self,pl,mode='normal'):
        """mode defines caching of power spectra for biases
        mode can be 'store_fid', 'use_fid' or normal"""
        bg=camb.get_background(self.pars)
        nz = 100 #4620 #number of steps to use for the radial/redshift integration
        kmax=10  #kmax to use
        #For Limber result, want integration over \chi (comoving radial distance), from 0 to chi_*.
        #so get background results to find chistar, set up arrage in chi, and calculate corresponding redshifts
        chistar = bg.conformal_time(0)- model.tau_maxvis.value
        chis = np.linspace(0,chistar,nz)
        zs=bg.redshift_at_comoving_radial_distance(chis)
        #zs=zs[::-1]
        Da=interp1d(zs,bg.comoving_radial_distance(zs))  
        Hi=interp1d(zs,1./(bg.h_of_z(zs)))   
        if (mode=='store_fid'):
            self.Da_fid=Da
            self.Hi_fid=Hi
            self.cpk_cached, self.mu_cached=[],[]
        pkl=[]
        for i,z in enumerate(self.zvals):
            if (mode=='use_fid'):
                cpk=self.cpk_cached[i]
                mu=self.mu_cached[i]
            else:
                kperp_t=self.kperp/Da(z)*self.Da_fid(z) ## we are observing radians, so..
                kpar_t=self.kpar/Hi(z)*self.Hi_fid(z)
                kt=np.sqrt(kperp_t**2+kpar_t**2)
                mu=kpar_t/kt
                cpk=[self.PK.P(z,k) for k in kt.flatten()]
                cpk=np.array(cpk).reshape(kt.shape)
                M=np.zeros(kt.shape)
                A=np.zeros(kt.shape)
                for j in range(self.Nkmu2_row):
                    for k in range(self.Nkmu2_col):
                        m=pl.value('Mkmu2_'+str(i)+str(j)+str(k))
                        a=pl.value('Akmu2_'+str(i)+str(j)+str(k))
                        if (a==0) and (m==0):
                            continue
                        X=(kt**j)*(mu**(2*k))
                        M+=X*m
                        A+=X*a
                cpk=A+cpk*(1+M)
                if (mode=='store_fid'):
                    self.cpk_cached.append(cpk)
                    self.mu_cached.append(mu)
            
            #f=self.growth_f(z,pl)
            #bpk=cpk*(pl.value('b_delta_'+str(i))+pl.value('b_eta_'+str(i))*f*mu**2)**2

            #Termporarily use the growth rate f(z) computed in CLASS for consistency.
            f_class=np.array((0.557204576947,
                              0.621162012915,
                              0.689749440772,
                              0.760259191738,
                              0.828475706028,
                              0.889167683497,
                              0.93742162152,
                              0.970389709813,
                              0.988403796107))
            bpk=cpk*(pl.value('b_delta_'+str(i))+pl.value('b_eta_'+str(i))*f_class[i]*mu**2)**2

            pkl.append(bpk)            
        return pkl

    def growth_f(self,z,pl):
        """Linear growth rate f(z) calculated using section 2.1 of
        arXiv:1405.1452"""
        om0=pl.value('omegac')+pl.value('omegab')
        Ez=np.sqrt(1-om0+om0*(1+z)**3.) #This is defined as H(z)/H(0)
        omz=om0*((1+z)**3.)/(Ez**2.)
        gamma=0.545
        f=omz**gamma
        return f


    def ComputeCosmo(self,pl):
        self.pars.set_cosmology(tau = pl.value('tau'), 
                                omch2 = pl.value('omegac'), 
                                H0 = None,
                                cosmomc_theta = pl.value('theta'), 
                                num_massive_neutrinos=1,
                                mnu = pl.value('mnu'),       
                                nnu = pl.value('Neff'),
                                standard_neutrino_neff = pl.value('Neff'),
                                #massless_neutrinos = pl.value('Neff')-1,  
                                ombh2 = pl.value('omegab'))
        #self.pars.set_dark_energy()
        self.pars.InitPower.set_params(As = pl.value('As'),
                                       ns = pl.value('ns'))
        #print ("Calling CAMB compute...",end='')

        #Get matter power spectrum interpolation object
        self.PK=camb.get_matter_power_interpolator(self.pars)  
        #print ("done")
