import numpy as np
from matplotlib import gridspec
import  csv
import pandas as pd
import matplotlib.pyplot as plt
import math as m
import utm
import sys
import itertools
import astropy.coordinates as coord
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, concatenate
from astropy.time import Time
import healpy as hp
from scipy.optimize import curve_fit
import seaborn as sn
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
from collections import Counter
from scipy import special
from scipy.stats import poisson

def reader(filename):
    df = pd.read_hdf(filename)
    return df

"""
after first inspction of the keys

Index(['TantraLines', 'TantraHits', 'Trigger3N', 'TriggerT3', 'Nrun',
       'TriggerCounter', 'FrameIndex', 'MCNGenEvent', 'MCNuType',
       'MCLeadingParticle', 'MCMultiplicity', 'Mestimator', 'TantraZenith',
       'TantraAzimuth', 'TantraAngularEstimator', 'TantraX', 'TantraY',
       'TantraZ', 'NOnTime', 'GridQuality', 'Lambda', 'AAZenith', 'AAAzimuth',
       'Beta', 'TrackLength', 'WeightAtmo', 'WeightAstro', 'TantraEnergy',
       'RunDurationYear', 'DateMJD', 'MCE', 'MCZenith', 'MCAzimuth', 'MCX',
       'MCY', 'MCZ', 'MCRho', 'w2', 'w3', 'MCRa', 'MCDec',
       'NEW_LIKELIHOOD_3D_ATMO', 'NEW_LIKELIHOOD_3D_ASTRO', 'IntegralCharge',
       'MeanCharge', 'StdCharge', 'BDT_default', 'BDT__cuts_1e2', 'TantraRho'],
      dtype='object')
"""
x=np.linspace(10**4,10**6,100) 
def limit_creation(E,norm,index):
    y=np.array(a*E**(-index+2))
    return y
    
if __name__ == "__main__":
    index_vector=np.array([2.0,2.2,2.4,2.5,2.7])
    norm_vector=np.array([2.4,2.3,2.0,1.9,1.8])#**10^-17 @ 100 TeV
    fig=plt.figure()
    ax=fig.add_subplot(111)
    for a in range(len(index_vector)):
        ax.plot(x,limit_creation(x,norm_vector[a]*10**-17,index_vector[a]))
        ax.set_yscale('log')
        ax.set_xscale('log')
    plt.show()
