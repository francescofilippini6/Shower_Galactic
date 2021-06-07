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
from mpl_toolkits.mplot3d import Axes3D
import threading

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
def cut_dataframe_bdts(df):
    #bdt_cut=0.33
    bdt_cut=0.12
    selectedDF=df[df['BDT__cuts_1e2'] < bdt_cut]
    return selectedDF

def type_appender(dfs,name):
    interaction_type=name.split('.')[0]
    print(interaction_type)
    y=[]
    for na in range(len(dfs['TantraLines'])):
        y.append(interaction_type)
    dfs['interaction_type']=y
    return dfs


def label_appender(dfaa,name):
    # dfaa = df[['TantraLines', 'TantraHits', 'Mestimator', 'TantraZenith',
    #           'TantraAzimuth', 'TantraAngularEstimator', 'TantraX', 'TantraY',
    #           'TantraZ','Lambda','Beta', 'TrackLength','TantraEnergy','TantraRho','IntegralCharge','MeanCharge', 'StdCharge',
    #            'GridQuality','AAZenith', 'AAAzimuth','Trigger3N', 'TriggerT3','NOnTime','NEW_LIKELIHOOD_3D_ATMO','MCE']]
    y=[]
    if 'numuCC' in name:
        y=np.ones(len(dfaa['TantraLines']))
        
    elif 'MUON' in name:
        y=np.ones(len(dfaa['TantraLines']))
    else:
        y=np.zeros(len(dfaa['TantraLines']))
    dfaa['label']=y
    return dfaa

if __name__ == "__main__":
    filename=sys.argv[1].split(',')
    listofdataframe=[]
    for counter,a in enumerate(filename):
        #keep an eye on the cut and on >< in the value selection
        print(a)
        df=reader(a)
        if "MUON" in a:
            print("Before cut",len(df['TantraX']))
            df1=label_appender(df,a)
            df1=type_appender(df1,a)
            listofdataframe.append(df1)
            print("after cut",len(df1['TantraX']))
        else:
            print("Before cut",len(df['TantraX']))
            df1=cut_dataframe_bdts(df)
            df1=label_appender(df1,a)
            df1=type_appender(df1,a)
            listofdataframe.append(df1)
            print("after cut",len(df1['TantraX']))
            
        
    result = pd.concat(listofdataframe)
    print(result.keys())
    print("final length",len(result['TantraX']))
    result.to_hdf('Prediction_MC_MUON_BDT_less0.12.h5', key='df', mode='w')
    
