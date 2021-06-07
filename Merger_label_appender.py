import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math as m
import sys
import utm
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
    dfAA = pd.read_hdf(filename)
    return dfAA

def cutter_interaction(dfcut,name):
    dfcut=dfcut[dfcut['interaction_type']==name]
    return dfcut

def selection(dfname):
   dfname=dfname[['Mestimator', 'TantraZenith', 'TantraAzimuth', 'Lambda', 'Beta',
       'WeightAtmo', 'TantraEnergy', 'RunDurationYear', 'DateMJD', 'MCE',
       'MCZenith', 'MCAzimuth', 'MCX', 'MCY', 'MCZ', 'MCRho', 'w2', 'w3',
                  'MCRa', 'MCDec', 'BDT_default', 'BDT__cuts_1e2', 'TantraRho','label']]
   return dfname

if __name__ == "__main__":
    filename=sys.argv[1]
    dfa=reader(filename)
    dfa=dfa.drop(['cont_label_pred','predicted_label'],axis=1)
    print(dfa.keys())
    originalfile=sys.argv[2]
    df1=reader(originalfile)
    print(df1.keys())
    interactions=np.unique(dfa['interaction_type'])
    listofdataframe=[]
    for inter in interactions:
        print(inter)
        dfb=cutter_interaction(dfa,inter)
        df2=cutter_interaction(df1,inter)
        print("length of intersection",len(dfb.index.intersection(df2.index)))
        comparison=np.array(dfb.index)==np.array(df2.index)
        print("Comparison, second check",comparison.all())
        dfsliced=df2.loc[np.array(dfb.index)]
        #dfb=selection(dfb)
        #dfsliced=selection(dfsliced)
        #print("Final comparison:",dfb.equals(dfsliced))
        dfb['predicted_dropout']=dfsliced['predicted_dropout']
        listofdataframe.append(dfb)
    dffinal= pd.concat(listofdataframe)
    print(dffinal.keys())
    #dffinal=dffinal.drop(['cont_label_pred','predicted_label'],axis=1)
    print(dffinal.keys())
    print("Writing to file")
    #dffinal.to_hdf('NewPrediction_NOTG_BDT_more_0.12.h5', key='df', mode='w')
    dffinal.to_hdf('Merged_NOTG_BDT_more_012.hdf', '/df')
        
    
