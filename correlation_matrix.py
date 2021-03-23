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
def correlation_matrix(df):
    #df=df.head(100)
    df1 = df[['TantraLines', 'TantraHits', 'Mestimator', 'TantraZenith',
              'TantraAzimuth', 'TantraAngularEstimator', 'TantraX', 'TantraY',
              'TantraZ','GridQuality', 'Lambda', 'AAZenith', 'AAAzimuth',
              'Beta', 'TrackLength','TantraEnergy','TantraRho','Trigger3N', 'TriggerT3','TriggerCounter','NOnTime']]
    corrMatrix = df1.corr()
    sn.heatmap(corrMatrix)
    plt.show()

def coordinate_plooter(df):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(df['TantraX'],df['TantraY'],df['TantraZ'], zdir='z', s=20, c=None, depthshade=True)
    plt.show()
    return fig

if __name__ == "__main__":
    filename=sys.argv[1]
    df=reader(filename)
    #coordinate_plooter(df)
    correlation_matrix(df)
    
