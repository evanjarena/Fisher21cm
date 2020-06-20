#
# Matrix for a tracer measuring Pk.
#

from __future__ import division, print_function
import sys
import numpy as np
import scipy.linalg as la
from scipy.interpolate import interp1d
from FishMat import FishMat
from ParameterVec import DefaultParamList, Parameter
from ClassWrap import PkDiffer_Class
from CambWrap import PkDiffer_Camb
import matplotlib.pyplot as plt

class TracerPk(FishMat):

    def __init__ (self, boltzmann_code, zvals, zmin, zmax, exp_name='None', SNR='None', kmax=0.5, dk=0.01, fsky=0.5, Nkmu2_row=3, Nkmu2_col=3):
        pl=DefaultParamList()
        self.N=len(pl)
        #ignorelist=['tau','As']
        self.ignorelist=[]

        self.boltzmann_code=boltzmann_code

        self.kvals=np.arange(dk,kmax+dk,dk)  
        self.Nk=len(self.kvals)
        self.dk=dk
        self.fsky=fsky
        self.kpar=np.outer(self.kvals,np.ones(self.Nk))
        self.kperp=self.kpar.T
        #plt.imshow(self.kperp)
        #plt.show()
        self.kt=np.sqrt(self.kpar**2+self.kperp**2)
        self.mu=self.kpar/self.kt
        self.Nkmu2_row=Nkmu2_row
        self.Nkmu2_col=Nkmu2_col

        # Find the redshift values to be used in normal mode calculation, calcNModes()
        zmax=[]
        zmin=[]
        for i in range(len(zvals)):
            if i == 0:
                zlow=0.
                zhigh=2.*zvals[i]
            else:
                zlow=zmax[i-1]
                zhigh=zvals[i]+zlow-zmin[i-1]
            zmin.append(zlow)
            zmax.append(zhigh)
        self.zvals=zvals
        self.zmax=zmax
        self.zmin=zmin
        assert(zmax[1] == zmin[2])
        self.Nz=len(zvals)

        self.exp_name=exp_name

        self.SNR=SNR[:,:self.Nk,:self.Nk]

        # Add bias parameters to parameter list
        for i,z in enumerate(self.zvals):
            pl.append(Parameter('b_delta_'+str(i), self.bias(z), ''))
            pl.append(Parameter('b_eta_'+str(i), self.biaseta(z), ''))
        self.Nwb=len(pl) # with biases

        # Add (a_zij)(kt^i)(mu^2)^j parameters to parameter list
        # Note that the fiducial a_zij values are zero
        for i in range(self.Nz):
            for j in range(self.Nkmu2_row):
                for k in range(self.Nkmu2_col):
                    pl.append(Parameter('Mkmu2_'+str(i)+str(j)+str(k),0,''))
                    pl.append(Parameter('Akmu2_'+str(i)+str(j)+str(k),0,''))
        self.Nwbkmu2=len(pl) # with biases and kmu2 parameters

        self.pl=pl

    #def Boltzmann_code(self):
    #    return self.boltzman_code

    def Pnoise(self,kpar,kperp,z):
        return 0

    def kNL(self,z):
        """Returns the comoving wavenumber on the threshold
        of the linear and non-linear regimes.
        """
        return min(0.5,0.04+0.016*(1+z)**2.2)

    def bias(self,z): 
        return 1.
        
    def biaseta(self,z):
        return 1.

    def kmu2(self):
        """Returns a Nkmu2_row x Nkmu2_col matrix of (kt^i)(mu^2)^j values.
        """
        kmu2=[]
        for i in range(self.Nkmu2_row):
            kmu2_row=[]
            for j in range(self.Nkmu2_col):
                kmu2_element=(self.kt**i)*(self.mu**2.)**j
                kmu2_row.append(kmu2_element)
            kmu2.append(kmu2_row)
        return kmu2

    def getInverseErrors(self):
        """Returns the inverse errors to be used in the Fisher matrix calculation.
        If SNR as a function of k_par and k_perp is specified for a given experiment,
        PkEI=1/\Delta P=SNR/P is returned.  
        If a specific experiment is not specified and hence no SNR is given,
        PkEI is calculated using the normal modes from calcNModes().
        """
        if (not hasattr(self,"nmodes")):
            self.calcNModes()
        PkC=self.PkDiffer.cube0
        PkEI=[] 
        if self.SNR=='None':
            for z,P,nm, in zip(self.zvals,PkC,self.nmodes):
                Pn=self.Pnoise(self.kpar,self.kperp,z)
                PkE=(P+Pn)**2/nm
                knl=self.kNL(z)
                PkE[np.where(self.kt>knl)]=1e30
                PkEI.append(1/PkE)
        else:
            for z,P,snr in zip(self.zvals, PkC, self.SNR):
                knl=self.kNL(z)
                PkE=P/snr
                PkE[np.where(self.kt>knl)]=1e30
                PkEI.append(1/PkE)
        return PkEI

    def calcNModes(self):
        self.nmodes=[]
        Da=self.PkDiffer.Da_fid
        for zlow,z,zhigh in zip(self.zmin,self.zvals,self.zmax):
            V=self.fsky*4*np.pi/3*(Da(zhigh)**3-Da(zlow)**3)
            Vk=2*np.pi*self.kperp*self.dk*self.dk
            cnm=V*Vk/(2*(2*np.pi)**3)
            self.nmodes.append(cnm)

    def calcFisher(self):
        """Calculates the Fisher matrix
        """
        print("Setting up Boltzmann code...") 

        if self.boltzmann_code=='Class':
            PkDiffer=PkDiffer_Class(self.pl,self.zvals, self.kvals, self.kperp, self.kpar, self.Nkmu2_row, self.Nkmu2_col)
        elif self.boltzmann_code=='Camb':
            PkDiffer=PkDiffer_Camb(self.pl,self.zvals, self.kvals, self.kperp, self.kpar, self.Nkmu2_row, self.Nkmu2_col)
        else:            
            print('WTF?')
        print('Your Boltzmann code of choice is ', self.boltzmann_code, '...')
        self.PkDiffer=PkDiffer

        PkEI=self.getInverseErrors()
        self.PkEI=PkEI
        eps=0.01
        Pderivs=[]
        print("Calculating derivatives... ", end='')
        for i1,p in enumerate(self.pl):
            print (" %s"%p.name, end='')
            sys.stdout.flush()
            if p.name not in self.ignorelist:
                Ders=self.PkDiffer.getDerivative(p,eps)
            else:
                Ders=None
            Pderivs.append(Ders)
        self.Pderivs=Pderivs
        print("")
        F1=np.zeros((self.Nwbkmu2,self.Nwbkmu2))
        print ("Getting fisher: ",end='')
        for i1,D1 in enumerate(Pderivs):
            print (" %i"%i1,end='')
            sys.stdout.flush()
            if D1 is None:
                continue
            for i2,D2 in enumerate(Pderivs):
                if (i2<i1):
                    continue
                if D2 is None:
                    continue
                for zi,z in enumerate(self.zvals):
                    v=(D1[zi]*PkEI[zi]*D2[zi]).sum() 
                    F1[i1,i2]+=v
                F1[i2,i1]=F1[i1,i2]
        F1+=np.diag([1e-30]*self.Nwbkmu2)

        # Fisher matrix with bias and kmu2 parametetrs
        F1_bkmu2=F1
        C_bkmu2=la.inv(F1_bkmu2)[:self.N,:self.N]
        F_bkmu2=la.inv(C_bkmu2)
        FishMat.__init__(self, self.pl[:self.N], F_bkmu2)
        if self.exp_name != 'None':
            FishMat.saveF(self, F_bkmu2, self.exp_name+'_bkmu2')

        # Fisher matrix with bias parameters only
        F1_b=F1[:self.Nwb,:self.Nwb]
        C_b=la.inv(F1_b)[:self.N,:self.N]
        F_b=la.inv(C_b)
        FishMat.__init__(self, self.pl[:self.N], F_b)
        if self.exp_name != 'None':
            FishMat.saveF(self, F_b, self.exp_name+'_b')

        # Print Fisher matrix with bias and kmu2 parameters
        print ('\n')
        print (F_bkmu2[:self.N,:self.N])


#        plt.figure()
#        plt.imshow(np.log(F1),interpolation='nearest')
#        plt.colorbar()
#        plt.show()


