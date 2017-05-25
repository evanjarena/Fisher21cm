from classy import Class
from parameters import Params
import numpy as np 
import h5py
from scipy.interpolate import interp1d

def power(paramvec, zs, ks):
    """Calculates the matter power spectrum P(k)
    """
        
    assert(len(paramvec)==8+len(zs))

    # Redshift bins
    zstr=",".join(map(str,zs))
    for z in zs:
        zstr+="%f"%z

    # 8 cosmological parameters
    pars = {
        'output': 'mPk',
        'P_k_max_h/Mpc': 20.,
        'tau_reio': paramvec[0],       
        'omega_cdm': paramvec[1],     
        'A_s': paramvec[2],
        '100*theta_s' : paramvec[3],           
        'N_ur': paramvec[4],  
        'N_ncdm': +1.,           
        'm_ncdm': paramvec[5],       
        'omega_b': paramvec[6],     
        'n_s': paramvec[7],          
        'z_pk' : zstr
    }

    # Bias parameters
    biases=paramvec[8:]

    # Create an instance of the CLASS wrapper
    cosmo = Class()

    # Set the parameters to the cosmological code
    cosmo.set(pars)

    # Run the whole code. Depending on your output, it will call the
    # CLASS modules more or less fast. For instance, without any
    # output asked, CLASS will only compute background quantities,
    # thus running almost instantaneously.
    # This is equivalent to the beginning of the `main` routine of CLASS,
    # with all the struct_init() methods called.
    cosmo.compute()

    #bg=cosmo.get_background()
    #bg_zs=bg['z']
    #Da=interp1d(bg_zs,bg['comov. dist.'])## cosmo.pk is actually all Mpc units
    #Hi=interp1d(bg_zs,1./(bg['H [1/Mpc]'])) # 

    # Growth rate
    # Use aproximate (non-interpolated) growth rate for now
    f=np.array((0.99279,
                0.977573,
                0.948137,
                0.907666,
                0.863469,
                0.791581,
                0.720339,
                0.675761,
                0.62456))

    # Compute the power spectrum  
    Pk=[]
    for zi, z in enumerate(zs):
        kpar=np.outer(np.ones(Nk),ks)
        kperp=np.outer(ks,np.ones(Nk))
        k=np.sqrt(kpar**2.+kperp**2.)

        mu=kpar/k
        beta=f[zi]/biases[zi]
    
        cpk=[cosmo.pk(k_val,z) for k_val in k.flatten()]
        cpk=np.array(cpk).reshape(k.shape)
        cpk*=(biases[zi]*(1+beta*mu**2.))**2.

        kmax=0.04+0.016*(1+z)**2.2
        cpk[np.where(k>kmax)]=1e30

        Pk.append(cpk)        
    return np.array(Pk)
