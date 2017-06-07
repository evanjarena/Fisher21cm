from classy import Class
from parameters import Params
from power_spectrum import power
import numpy as np 
import h5py
import pickle

# Data cube from Cosmic Visions simulation
da=h5py.File("sn_lowz_expA_50K.h5")

# Redshift bins
zs=[]
for i in range(9):
    czs=da['sn_band_'+str(i)].attrs['z']
    zs.append(czs)
Nz=len(zs)

# Signal to noise ratio (SNR) data
snr=[]
for i in range(9):
    csnr=da['sn_band_'+str(i)].value
    snr.append(csnr)
snr=np.array(snr)

# Comoving wave numbers
ks=da['k_bin'].value
# The length of the kbin in the hp5y file is 501  
# We need to make it 500 to match the length of the SNR data
#Let us start at 0.01 and end at 5.00
ks=ks[1:]
Nk=len(ks)

# Biases
# Use approximate (non-interpolated) biases for now
biases=np.array((3.7753,
                 2.4726,
                 1.98777,
                 1.7671,
                 1.6526,
                 1.56537,
                 1.50723,
                 1.47815,
                 1.44907))

# Call instance of Params class
P=Params(biases)
Np=len(P.default_params())

# Define power spectrum
defpower=power(P.default_params(),zs,ks)

# Calculate derivatives
ders=[]
steps=P.default_params()*0.01
for i in range(Np):
    step=steps[i]
    params_plus=P.perturbed_params(i,step)
    params_minus=P.perturbed_params(i,-step)
    #powerplus=np.log(power(params_plus,zs,ks))
    #powerminus=np.log(power(params_minus,zs,ks))
    powerplus=power(params_plus,zs,ks)
    powerminus=power(params_minus,zs,ks)
    deriv=(powerplus-powerminus)/(2*step)/defpower
    #print 'parameter', i, 'derivatives are', deriv
    ders.append(deriv)

# Now find the Fisher and covariance matrices
Ft=np.zeros((Np,Np))
for zi in range(Nz):
    Fc=np.zeros((Np,Np))
    for i in range(Np):
        for j in range(i,Np):
            Fc[i,j]=(ders[i][zi,:,:]*ders[j][zi,:,:]*snr[zi]**2.).sum()
            Fc[j,i]=Fc[i,j]
    Ft+=Fc

# Save Fisher matrix to a .txt file for quick glance
# Let delimieter be & in order to more easily put matrix in LaTeX Tabular environment
# (Add full conversion to LaTeX table)
np.savetxt('fisher21cm_matrix.txt', Ft, delimiter=' & ', fmt='%0.4f')

# Save Fisher matrix into a pickle file to load into other modules
pickle.dump(Ft, open('fisher21cm_matrix.p', 'wb'))
