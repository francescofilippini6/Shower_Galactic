import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math as m
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

def cut_on_prediction(df):
    df=df[df['predicted_label']<0.1]
    return df

def weight_astro_spectrum(df):
    spectral_index = 2.0
    spectral_norm = 6 * 10**-12 * 10**(3*spectral_index)  #Neutrino Energy in GeV
    #new_w3_weight=np.array(df['w2'])*spectral_norm*np.power(np.array(np.power(10,df['MCE'])),-spectral_index)
    my_spectral_index = 2.5
    my_spectral_norm = 3.5 * 10**-12 * 10**(3*my_spectral_index)  #Neutrino Energy in GeV
    #imporvements of cpu time of a factor 1000 minimum !!!!
    new_w3_weight=np.array(df['w2'])*my_spectral_norm*np.power(np.array(np.power(10,df['MCE'])),-my_spectral_index)
    #error_norm=(4.9 *10**-12 * 10**(3*my_spectral_index))*np.power(np.array(np.power(10,df['MCE'])),-my_spectral_index)*np.array(df['w2']) #energy of nu  in GeV
    df['new_w3'] = np.array(new_w3_weight)
    #df['error_norm'] = np.array(error_norm)
    print("done")
    return df


def plot_reco_energy_distro(df):
    livetime=3012/365
    fig = plt.figure()
    ax=fig.add_subplot(111)
    #ax.hist(df['TantraEnergy'],bins=50,histtype='step',weights=np.array(df['new_w3'])*livetime,label='nu cosmic total')
    print("first cut")
    df=df[df['predicted_label']<0.1]
    print("first histo")
    ax.hist(df['TantraEnergy'],bins=50,histtype='step',weights=np.array(df['WeightAtmo'])*livetime,label='nu atmo, prediction 0')
    print("second histo")
    ax.hist(df['TantraEnergy'],bins=50,histtype='step',weights=np.array(df['new_w3'])*livetime,label='nu cosmic, prediction 0, -2.4')
    print("second cut")
    df1=df[df['predicted_label']>0.8]
    print("third histo")
    ax.hist(df1['TantraEnergy'],bins=50,histtype='step',weights=np.array(df1['WeightAtmo'])*livetime,label='nu atmo, prediction 1')
    ax.set_xlabel('log10(E/GeV)')
    ax.set_ylabel(r'$\frac{dN}{dE}$')
    ax.set_yscale('log')
    ax.legend()
    plt.show()
    return fig

if __name__ == "__main__":
    filename=sys.argv[1]
    df=reader(filename)
    #print("LENGTH:",len(df['TriggerT3']))
    #print(Counter(df['label']))
    #print(Counter(df['interaction_type']))
    print(df.keys())
    #df1=weight_astro_spectrum(cut_on_prediction(df))
    df1=weight_astro_spectrum(df)
    print("plotting")
    plot_reco_energy_distro(df1)
