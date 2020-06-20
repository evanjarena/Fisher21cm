#
# Super simple Fisher matrix class.
#

from __future__ import division, print_function
import numpy as np
import scipy.linalg as la
from ParameterVec import DefaultParamList
import pickle

class FishMat(object):

    def __init__ (self, plist, mat):
        self.plist=plist
        self.F=mat
        self.N=len(self.F)
        if (self.N!=len(plist)):
            print ("Bad FishMat constructor")
            stop()
        self.calcC()

    def calcC(self):
        self.C=la.inv(self.F+np.diag([1e-30]*self.N))

    def saveF(self, F, exp_name):
        """Save Fisher matrix into a pickle file to load into other modules
        """
        pickle.dump(F, open('FishData/'+str(exp_name)+'_fisher_matrix.p', 'wb'))
    
    def addF(self,F):
        if (type(F)==list):
            Fl=F
        else:
            Fl=[F]
        for Fa in Fl:
            if self.plist.nameList()==Fa.plist.nameList():
                #print(self.F)
                #print(Fa)
                self.F+=Fa.F
            else:
                print(self.plist.nameList())
                print(Fa.plist.nameList())
                print("Params don't match")
                stop()
        self.calcC()
        
    def error(self,pname):
        i=self.plist.nameList().index(pname)
        return np.sqrt(self.C[i,i])
    
    def marginalise(self,pname_list):
        found=0
        outndx=[]
        for i,n in self.plist:
            if n in pname_list:
                found+=1
            else:
                outndx.append(i)
        if found!=len(pname_list):
            print("Didn't find all parameters in marginalise")
            stop()
        Nn=len(outndx)
        outC=np.zeros((Nn,Nn))
        for inew,iold in enumerate(outndx):
            outC[inew,:]=self.C[iold,outndx]
            outC[:,inew]=self.C[outndx,iold]
        self.C=outC
        self.F=la.inv(outC)
        
