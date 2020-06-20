#
# Joel's Fisher matrix for CMBS4.
#

from __future__ import division, print_function
import numpy as np
import scipy.linalg as la
from FishMat import FishMat
from ParameterVec import DefaultParamList

class S4Fisher(FishMat):
    
    def __init__ (self):
        m=np.loadtxt('FishData/s4fisher05042017.txt')
        pl=DefaultParamList()
        assert(pl.nameList()==['tau','omegac','As','theta','Neff','mnu','omegab','ns'])
        FishMat.__init__(self,pl,m)
        FishMat.saveF(self, m, 'CMBS4')

S4Fisher()
