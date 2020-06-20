#
# Create triangle plots and deal with multiple Fisher matrices.
#

#from __future__ import division, print_function
from ParameterVec import DefaultParamList
import numpy as np
import pickle
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.patches as mpatches

class FisherAnalysis(object):

    def __init__(self, params, param_names):
        self.params=params
        self.Np=len(params[0:8]) # 8 parameter model
        self.param_names=param_names
        self.fisher_matrix=0.
        self.covariance_matrix=0.
        self.fms = []
        self.Nfms = len(self.fms)
        self.cms=[]
        self.Ncms=len(self.cms)
        self.clrs = ["cornflowerblue", "g", "r", "g", "r", "r"]
        self.linestyles = ['solid', 'solid', 'solid', 'solid']
        self.linewidths = [1.5, 2.0, 2.5, 1.5]
        self.fills = [True, True, True, True]
        self.alphas = [0.75, 0.35, 1.0, 0.35]

    def import_fisher(self, filename):
        """ Import a saved fisher matrix using pickle.
        Then compute the covariance matrix, marginalize over any bias
        parameters in the matrix, and return the NpxNp Fisher matrix.
        """
        self.fisher_matrix=pickle.load(open(filename, 'rb'))
        cov=np.linalg.inv(self.fisher_matrix)
        cov=cov[0:self.Np,0:self.Np]
        self.fisher_matrix=np.linalg.inv(cov)
        return self.fisher_matrix

    def covariance(self, fisher_matrix):
        """Compute the covariance matrix
        """
        self.covariance_matrix=np.linalg.inv(fisher_matrix)
        return self.covariance_matrix

    def covariance_array(self, *covariance_matrices):
        """Creates an array of covariance matrices for combined visualization
        """
        self.cms=covariance_matrices
        self.Ncms=len(covariance_matrices)
        return self.cms

    def error(self, covariance_matrix):
        """Computes the marginalized 1-sigma uncertainty on each parameter
        """
        errors=[]
        for i in range(self.Np):
            err=np.sqrt(covariance_matrix[i,i])
            errors.append(err)
        return errors
        

    def marginalize(self, covariance_matrix, i, j):
        """Compute and return a new covariance matrix after marginalizing over all
        other parameters not in the list param_list
        """
        param_list=[i,j]

        mlist=[]
        for k in range(self.Np):
            if k not in param_list:
                mlist.append(k)
        
        trun_cov_matrix=covariance_matrix        
        trun_cov_matrix=np.delete(trun_cov_matrix, mlist, 0)
        trun_cov_matrix=np.delete(trun_cov_matrix, mlist, 1)
        return trun_cov_matrix

    def error_ellipse(self, covariance_matrix, i, j, nstd=1, space_factor=5., clr='b', alpha=0.5, lw=1.5):
        """Compute 2D marginalized confidence ellipses
        """
        def eigsorted(cov):
            vals, vecs=np.linalg.eigh(cov)
            order=vals.argsort()[::1]
            return vals[order], vecs[:,order]

        marginal_covariance_matrix=self.marginalize(covariance_matrix, i,j)
        vals, vecs=eigsorted(marginal_covariance_matrix)
        theta=np.degrees(np.arctan2(*vecs[:,0][::-1]))
        width, height=2*nstd*np.sqrt(np.absolute(vals))
        xypos=[self.params[i], self.params[j]]
        ellip=Ellipse(xy=xypos, width=width, height=height, angle=theta, color=clr, alpha=alpha)
        ellip.set_linewidth(lw)
        ellip_vertices=ellip.get_verts()
        xl=[ellip_vertices[k][0] for k in range(len(ellip_vertices))]
        yl=[ellip_vertices[k][1] for k in range(len(ellip_vertices))]
        dx=(max(xl)-min(xl))/space_factor
        dy=(max(yl)-min(yl))/space_factor
        xyaxes=[min(xl)-dx, max(xl)+dx, min(yl)-dy, max(yl)+dy]

        return ellip, xyaxes

    def oneD_constraints(self, covariance_matrix, i, j, clr='b', lw=1.5):
        """Compute 1D marginalized constraints
        """
        # Ensure that we are computing the diagonal of the triangle plot
        assert(i==j)
         
        marginal_covariance_matrix=self.marginalize(covariance_matrix, i, j)
        xpos=self.params[i]

        sig=np.sqrt(np.absolute(marginal_covariance_matrix))
        sig=sig.item()
        xx=np.linspace(xpos-20.*sig, xpos+20.*sig, 4000)
        yy=1./np.sqrt(2.*np.pi*sig**2.) * np.exp(-0.5 * ((xx-xpos)/sig)**2.)
        yy/=np.max(yy)
        
        return yy, xx

    def plot_error_ellipse(self, covariance_matrix, i, j, xyaxes_input=0, nstd=1, clr='b', alpha=0.5, lw=1.5):
        """Plot the 2D marginalized error ellipses
        """
        ax = plt.gca()
        errorellipse, xyaxes=self.error_ellipse(covariance_matrix, i, j, nstd=nstd, clr=clr, alpha=alpha, lw=lw)
        ax.add_artist(errorellipse)
        if (xyaxes_input!=0):
            ax.axis(xyaxes_input)
        else:
            ax.axis(xyaxes)

    def plot_oneD_constraints(self, covariance_matrix, i, j, clr='b', lw=1.5):
        """Plot the 1D marginalized Gaussians
        """
        ax=plt.subplot(111)
        y,x=self.oneD_constraints(covariance_matrix, i, j, clr=clr, lw=lw)
        ax.plot(x,y)
        plt.show()

    def plot_error_matrix(self, covariance_matrix, filename, nstd=1, nbinsx=6, nbinsy=6):
        """Create a triangle plot of 2D error ellipses and 1D Gaussians
        given the list of parameters provided
        """
        Np=self.Np

        f, allaxes = plt.subplots(Np, Np, sharex="col", figsize=(10,10))#, sharey="row")

        for j in range(Np):
            for i in range(Np):
                # Off-diagonal 2D marginalized plots
                if (j>i):
                    # 1-sigma ellipse
                    errorellipse, xyaxes=self.error_ellipse(covariance_matrix, i,j, nstd=nstd, clr="cornflowerblue", alpha=0.75)
                    # 2-sigma ellipse
                    ere2, xyaxes2 = self.error_ellipse(covariance_matrix, i, j, nstd=2*nstd, clr="cornflowerblue", alpha=0.35)
                    
                    # Define axes
                    jp=i
                    ip=j
                    axis=allaxes[ip][jp]
                    axis.locator_params(axis='x', nbins=nbinsx)
                    axis.locator_params(axis='y', nbins=nbinsy)
                    axis.add_artist(ere2)
                    axis.add_artist(errorellipse)
                    axis.axis(xyaxes2)
                    if (i==0):
                        axis.set_ylabel(self.param_names[j], fontsize=14)
                    if (i>=0):
                        axis.set_xlabel(self.param_names[i], fontsize=14)
                    
                    # Make ticks visible and make sure they are not too crowded
                    axis.ticklabel_format(style='sci', axis='both', scilimits=(-3,3), useOffset=False)
                    for label in axis.get_xticklabels():
                        label.set_rotation(-90) 
                    axis.tick_params(labelsize=8)                       # Tick values
                    axis.yaxis.get_offset_text().set_fontsize(8)        # Scientific notation floats
                    axis.xaxis.get_offset_text().set_fontsize(8)        # Scientific notation floats
                    if i != 0:
                        for tick in axis.yaxis.get_major_ticks():
                            tick.label1.set_visible(False)
                """            
                if (j<i):
                    # 1-sigma ellipse
                    errorellipse, xyaxes=self.error_ellipse(covariance_matrix, i,j, nstd=nstd, clr="cornflowerblue", alpha=.3)
                    # 2-sigma ellipse
                    ere2, xyaxes2 = self.error_ellipse(covariance_matrix, i, j, nstd=2*nstd, clr="cornflowerblue", alpha=.75)
                    
                    # Define axes
                    jp=i
                    ip=j
                    axis=allaxes[ip][jp]
                    axis.locator_params(axis='x', nbins=nbinsx)
                    axis.locator_params(axis='y', nbins=nbinsy)
                    axis.add_artist(ere2)
                    axis.add_artist(errorellipse)
                    axis.axis(xyaxes2)
                    if (i==0):
                        axis.set_ylabel(self.param_names[j], fontsize=14)
                    if (i>=0):
                        axis.set_xlabel(self.param_names[i], fontsize=14)
                """
                # On-diagonal 1D marginalized plots
                if (j==i):
                    y, x=self.oneD_constraints(covariance_matrix, i,j, clr="cornflowerblue")
                    
                    jp=i
                    axis=allaxes[jp][jp]

                    axis.locator_params(axis='x', nbins=nbinsx)
                    axis.locator_params(axis='y', nbins=nbinsy)

                    if (i==0):
                        axis.set_ylabel(self.param_names[j], fontsize=14)
                    if (i>=0):
                        axis.set_xlabel(self.param_names[i], fontsize=14)

                    axis.ticklabel_format(style='sci', axis='both', scilimits=(-3,3), useOffset=False)
                    for label in axis.get_xticklabels():
                        label.set_rotation(-90) 
                    axis.tick_params(labelsize=8)                       # Tick values
                    axis.yaxis.get_offset_text().set_fontsize(8)        # Scientific notation floats
                    axis.xaxis.get_offset_text().set_fontsize(8)        # Scientific notation floats
                    if i != 0:
                        for tick in axis.yaxis.get_major_ticks():
                            tick.label1.set_visible(False)

                    axis.plot(x,y)

        # Hide empty plots above the diagonal
        for i in range(Np):
            for j in range(Np):
                if (j>i):
                    allaxes[i][j].axis('off')
                
        # Tight layout
        plt.tight_layout()

        # Remove whitespace between subplots
        plt.subplots_adjust(wspace=0, hspace=0)

        # Save figure
        plt.savefig(filename+'.png', format='png', dpi=1000) 
        plt.show()

    def plot_error_matrix_combined(self, covariance_matrix_array, filename, nstd=1, nbinsx=6, nbinsy=6):
        """Create a triangle plot for comparison of multiple experiments
        Some examples:
        This allows one to plot Experiment A in one color versus Exp. B in another color
        This also allows one to plot (Exp. A + Exp. B) versus Exp. A
        """
        Np=self.Np

        f, allaxes = plt.subplots(Np, Np, sharex="col", figsize=(10,10))#, sharey="row")

        for j in range(Np):
            for i in range(Np):
                if (j>i):
                    cmcnt=0
                    for c in range(len((covariance_matrix_array))):
                        errorellipse, xyaxes=self.error_ellipse(covariance_matrix_array[c], i,j, 
                                                                nstd=nstd, clr=self.clrs[cmcnt], alpha=self.alphas[cmcnt])
                        ere2, xyaxes2 = self.error_ellipse(covariance_matrix_array[c], i, j, 
                                                           nstd=2*nstd, clr=self.clrs[cmcnt], alpha=self.alphas[cmcnt]/2.)
                        errorellipse.set_linestyle(self.linestyles[cmcnt])
                        errorellipse.set_linewidth(self.linewidths[cmcnt])
                        errorellipse.set_fill(self.fills[cmcnt])
                        ere2.set_linestyle(self.linestyles[cmcnt])
                        ere2.set_linewidth(self.linewidths[cmcnt])
                        ere2.set_fill(self.fills[cmcnt])
                        
                        jp=i
                        ip=j
                        axis=allaxes[ip][jp]
                        
                        axis.add_artist(ere2)
                        axis.add_artist(errorellipse)
                        axis.axis(xyaxes2)
                        if cmcnt == 0:
                            axis.axis(xyaxes)
                        cmcnt=cmcnt+1
                        
                    axis.locator_params(axis='x', nbins=nbinsx)
                    axis.locator_params(axis='y', nbins=nbinsy)
                    if (i==0):
                        axis.set_ylabel(self.param_names[j], fontsize=14)
                    if (i>=0):
                        axis.set_xlabel(self.param_names[i], fontsize=14)

                    axis.ticklabel_format(style='sci', axis='both', scilimits=(-3,3), useOffset=False)
                    for label in axis.get_xticklabels():
                        label.set_rotation(-90) 
                    axis.tick_params(labelsize=8)                       # Tick values
                    axis.yaxis.get_offset_text().set_fontsize(8)        # Scientific notation floats
                    axis.xaxis.get_offset_text().set_fontsize(8)        # Scientific notation floats
                    if i != 0:
                        for tick in axis.yaxis.get_major_ticks():
                            tick.label1.set_visible(False)

                if (j==i):
                    cmcnt=0
                    for c in range(len((covariance_matrix_array))):
                        y, x=self.oneD_constraints(covariance_matrix_array[c], i,j, clr=self.clrs[cmcnt])
                        jp=i
                        axis=allaxes[jp][jp]
                        axis.plot(x,y)
                        cmcnt=cmcnt+1

                    axis.locator_params(axis='x', nbins=nbinsx)
                    axis.locator_params(axis='y', nbins=nbinsy)
                    if (i==0):
                        axis.set_ylabel(self.param_names[j], fontsize=14)
                    if (i>=0):
                        axis.set_xlabel(self.param_names[i], fontsize=14)

                    axis.ticklabel_format(style='sci', axis='both', scilimits=(-3,3), useOffset=False)
                    for label in axis.get_xticklabels():
                        label.set_rotation(-90) 
                    axis.tick_params(labelsize=8)                       # Tick values
                    axis.yaxis.get_offset_text().set_fontsize(8)        # Scientific notation floats
                    axis.xaxis.get_offset_text().set_fontsize(8)        # Scientific notation floats
                    if i != 0:
                        for tick in axis.yaxis.get_major_ticks():
                            tick.label1.set_visible(False)
                            
        for i in range(Np):
            for j in range(Np):
                if (j>i):
                    allaxes[i][j].axis('off')

        plt.tight_layout
        plt.subplots_adjust(wspace=0, hspace=0)

        #Save figure
        plt.savefig(filename+'.png', format='png', dpi=1000) 
        plt.show()

#---------------------------------------------------------------------------------------------------------

# Get parameter values and parameter names
pl=DefaultParamList()
P=pl.valueList()
N=pl.nameLaTeX()

##########################################################################################################
# Plots for Fisher matrices that include both biases and kmu2 matrix (both of which are marginalized over)
##########################################################################################################

# Get 21cm plot
F_21cm=FisherAnalysis(P,N).import_fisher('FishData/21cm_bkmu2_fisher_matrix.p')
C_21cm=FisherAnalysis(P,N).covariance(F_21cm)
FisherAnalysis(P,N).plot_error_matrix(C_21cm, 'FishData/21cm_bkmu2_triangle')

#Get CMBS4 plot
F_S4=FisherAnalysis(P,N).import_fisher('FishData/CMBS4_fisher_matrix.p')
C_S4=FisherAnalysis(P,N).covariance(F_S4)
FisherAnalysis(P,N).plot_error_matrix(C_S4, 'FishData/CMBS4_triangle')

##Plot of 21cm vs CMBS4
#C_21cm_vs_S4=FisherAnalysis(P,N).covariance_array(C_21cm, C_S4)
#FisherAnalysis(P,N).plot_error_matrix_combined(C_21cm_vs_S4, '21cm_vs_CMBS4_triangle')

##Plot of 21cm vs 21cm+CMBS4
#C_21cm_and_C_S4=FisherAnalysis(P,N).covariance(F_21cm+F_S4)
#C_21cm_vs_C_21cm_and_S4=FisherAnalysis(P,N).covariance_array(C_21cm, C_21cm_and_C_S4)
#FisherAnalysis(P,N).plot_error_matrix_combined(C_21cm_vs_C_21cm_and_S4, '21cm_vs_21cmCMBS4_triangle')

#Plot of CMBS4 vs 21cm+CMBS4
C_21cm_and_C_S4=FisherAnalysis(P,N).covariance(F_21cm+F_S4)
C_S4_vs_C_21cm_and_S4=FisherAnalysis(P,N).covariance_array(C_S4, C_21cm_and_C_S4)
FisherAnalysis(P,N).plot_error_matrix_combined(C_S4_vs_C_21cm_and_S4, 'FishData/CMBS4_vs_21cm_bkmu2_CMBS4_triangle')

# Get 1-sigma errors
print '\n'
print 'The following are the 1-sigma errors on parameters'
print 'tau, omegac, As, theta, Neff, mnu, omegab, ns'
print 'derived from matrices with both bias and kmu2 parameters'
print 'for the following experiment combinations:'
print '21cm errors are', FisherAnalysis(P,N).error(C_21cm)                   # 21cm alone
print 'CMBS4 errors are', FisherAnalysis(P,N).error(C_S4)                    # CMBS4 alone
print '21cm+CMBS4 errors are', FisherAnalysis(P,N).error(C_21cm_and_C_S4)    # 21cm + CMBS4 

##########################################################################################################
# Plots for Fisher matrices that include only biases (which are marginalized over) and not the kmu2 matrix 
##########################################################################################################

# Get 21cm plot
F_21cm=FisherAnalysis(P,N).import_fisher('FishData/21cm_b_fisher_matrix.p')
C_21cm=FisherAnalysis(P,N).covariance(F_21cm)
FisherAnalysis(P,N).plot_error_matrix(C_21cm, 'FishData/21cm_b_triangle')

#Get CMBS4 plot
F_S4=FisherAnalysis(P,N).import_fisher('FishData/CMBS4_fisher_matrix.p')
C_S4=FisherAnalysis(P,N).covariance(F_S4)
FisherAnalysis(P,N).plot_error_matrix(C_S4, 'FishData/CMBS4_triangle')

##Plot of 21cm vs CMBS4
#C_21cm_vs_S4=FisherAnalysis(P,N).covariance_array(C_21cm, C_S4)
#FisherAnalysis(P,N).plot_error_matrix_combined(C_21cm_vs_S4, '21cm_vs_CMBS4_triangle')

##Plot of 21cm vs 21cm+CMBS4
#C_21cm_and_C_S4=FisherAnalysis(P,N).covariance(F_21cm+F_S4)
#C_21cm_vs_C_21cm_and_S4=FisherAnalysis(P,N).covariance_array(C_21cm, C_21cm_and_C_S4)
#FisherAnalysis(P,N).plot_error_matrix_combined(C_21cm_vs_C_21cm_and_S4, '21cm_vs_21cmCMBS4_triangle')

#Plot of CMBS4 vs 21cm+CMBS4
C_21cm_and_C_S4=FisherAnalysis(P,N).covariance(F_21cm+F_S4)
C_S4_vs_C_21cm_and_S4=FisherAnalysis(P,N).covariance_array(C_S4, C_21cm_and_C_S4)
FisherAnalysis(P,N).plot_error_matrix_combined(C_S4_vs_C_21cm_and_S4, 'FishData/CMBS4_vs_21cm_b_CMBS4_triangle')

# Get 1-sigma errors
print '\n'
print 'The following are the 1-sigma errors on parameters'
print 'tau, omegac, As, theta, Neff, mnu, omegab, ns'
print 'derived from matrices with bias parameters only (no kmu2)'
print 'for the following experiment combinations:'
print '21cm errors are', FisherAnalysis(P,N).error(C_21cm)                   # 21cm alone
print 'CMBS4 errors are', FisherAnalysis(P,N).error(C_S4)                    # CMBS4 alone
print '21cm+CMBS4 errors are', FisherAnalysis(P,N).error(C_21cm_and_C_S4)    # 21cm + CMBS4
