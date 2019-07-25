###### Import stuff ###### 
import numpy as np
import astropy.units as unt
#import matplotlib.pyplot as plt
from astropy.cosmology import WMAP9 as cosmo
import sys
from scipy import interpolate
from scipy.integrate import nquad, quad
#%matplotlib inline
from scipy.constants import c, G
from astropy.constants import M_sun


def draw_powerlaw_samples(**kwargs):
    '''
    Yields random masses, with the first component drawn from
    the Salpeter initial mass function distribution and the
    second mass drawn uniformly between min_mass and the mass of
    the first component.
    '''
 
    nsamples = kwargs.get('nsamples', 1)
    min_mass = kwargs.get('min_mass', 5.)
    max_mass = kwargs.get('max_mass', 50.)
    max_mtotal = kwargs.get('max_mtotal', 2 * max_mass)
    alpha = kwargs.get('alpha', 2.35)    
        
    a = (max_mass/min_mass)**(-alpha + 1.0) - 1.0
    beta = 1.0 / (-alpha + 1.0)
    
    k = nsamples * int(2.0 + np.log(1 + 100./nsamples))
    aa = min_mass * (1.0 + a * np.random.random(k))**beta
    bb = np.random.uniform(min_mass, aa, k) 
    
    idx = np.where(aa + bb < max_mtotal)
    m1, m2 = (np.maximum(aa, bb))[idx], (np.minimum(aa, bb))[idx]
    
    m1.resize(nsamples)
    m2.resize(nsamples)    
    
    if m1[-1] == 0:
        sys.exit("Rejection sampling including zeros in the population masses!")    
        
    return m1, m2


def draw_flatlog_samples(**kwargs):
    '''
    Yields random masses drawn uniformly-in-log between min_mass
    and max_mass, discarding those with a total mass exceeding
    max_mtotal.
    '''
    #PDF doesnt match with sampler
    nsamples = kwargs.get('nsamples', 1)
    min_mass = kwargs.get('min_mass', 5.)
    max_mass = kwargs.get('max_mass', 50.)
    max_mtotal = kwargs.get('max_mtotal', 2 * max_mass)  
    flatlogmin = np.log(min_mass)
    flatlogmax = np.log(max_mass)
    
    k = nsamples * int(2.0 + np.log(1 + 100./nsamples))
    aa = np.exp(np.random.uniform(flatlogmin, flatlogmax, k))
    bb = np.exp(np.random.uniform(flatlogmin, flatlogmax, k))
  
    idx = np.where(aa + bb < max_mtotal)
    m1, m2 = np.maximum(aa[idx], bb[idx]), np.minimum(aa[idx], bb[idx])
    
    m1.resize(nsamples)
    m2.resize(nsamples)    
    
    if m1[-1] == 0:
        sys.exit("Rejection sampling including zeros in the population masses!") 
        
    return m1, m2



