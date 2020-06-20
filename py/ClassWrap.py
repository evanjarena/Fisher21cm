#
# Wrapper for Class.
#

from __future__ import division, print_function
from classy import Class
import copy
from scipy.interpolate import interp1d
from ParameterVec import DefaultParamList, ParamList, Parameter
import sys
import numpy as np

class PkDiffer_Class:

    def __init__ (self,pl, zvals, kvals, kperp, kpar, Nkmu2_row, Nkmu2_col):
        """ returns a list of Pks, each list containins 3D PS"""
        self.zstr=",".join(map(str,zvals+[zvals[-1]+2]))
        self.kvals=kvals
        self.kperp=kperp
        self.kpar=kpar
        self.zvals=zvals
        self.Nkmu2_row=Nkmu2_row
        self.Nkmu2_col=Nkmu2_col
        self.plist=copy.deepcopy(pl)
        self.cosmo = Class()
        self.ComputeCosmo(pl)
        bg=self.cosmo.get_background()
        zs=bg['z']
        zs=zs[::-1]
        Da=interp1d(zs,bg['comov. dist.'])  # cosmo.pk is actually all Mpc units
        Hi=interp1d(zs,1./(bg['H [1/Mpc]']))  
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
        #dp,dm=de
        #toret=(dp-dm)/(2*step)
        toret=[(dp-dm)/(2*step) for dp,dm in zip(de[0],de[1])]
        return toret
        
    def getCube(self,pl,mode='normal'):
        """mode defines caching of power spectra for biases
        mode can be 'store_fid', 'use_fid' or normal"""
        bg=self.cosmo.get_background()
        zs=bg['z']
        zs=zs[::-1]
        Da=interp1d(zs,bg['comov. dist.'])## cosmo.pk is actually all Mpc units
        Hi=interp1d(zs,1./(bg['H [1/Mpc]'])) # 
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
                #mu=self.kpar/np.sqrt(self.kperp**2+self.kpar**2)
                #[print(k) for k in kt.flatten()]
                #print(len(kt.flatten()))
                cpk=[self.cosmo.pk(k,z) for k in kt.flatten()]
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
            
            f=self.growth_f(z)
            bpk=cpk*(pl.value('b_delta_'+str(i))+pl.value('b_eta_'+str(i))*f*mu**2)**2
            pkl.append(bpk)
            
        return pkl

    def growth_f(self,z):
        da=0.01
        a=1./(1.+z)
        gp,g,gm=[self.cosmo.scale_independent_growth_factor(1./ia-1.) for ia in [a+da,a,a-da]]
        f=a*(gp-gm)/(2*g*da)
        return f

    def ComputeCosmo(self,pl):
        #del self.cosmo
        #self.cosmo = Class()
        pars = {
                'output': 'mPk',
                'P_k_max_h/Mpc': self.kvals[-1]+3.0,
                'tau_reio': pl.value('tau'),
                'omega_cdm': pl.value('omegac'),     
                'A_s': pl.value('As'),     
                '100*theta_s' : 100*pl.value('theta'),     
                'N_ur': pl.value('Neff')-1,  
                'N_ncdm': 1.0,           
                'm_ncdm': pl.value('mnu'),       
                'omega_b': pl.value('omegab'),     
                'n_s': pl.value('ns'),          
                'z_pk' : self.zstr,
        }
        self.cosmo.set(pars)
        #print ("Calling class compute...",end='')
        self.cosmo.compute()
        #print ("done")

        
 
