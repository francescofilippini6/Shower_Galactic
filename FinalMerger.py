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

   
if __name__ == "__main__":
    old_nomu=sys.argv[1]
    muon=sys.argv[2]
    listofdataframe=[]
    dfold=reader(old_nomu)
    df1=dfold.drop(['predicted_label'],axis=1)
    labels=np.array(df1['cont_label_pred'])
    df1=df1.drop(['cont_label_pred'],axis=1)
    df1['predicted_label']=labels
    df1=df1.drop(['Mestimator','MCX', 'MCY', 'MCZ', 'MCRho','MCRa', 'MCDec', 'BDT_default'],axis=1)
    listofdataframe.append(df1)
    df2=reader(muon)
    df2=df2.drop(['Mestimator','MCX', 'MCY', 'MCZ', 'MCRho','MCRa', 'MCDec', 'BDT_default'],axis=1)
    listofdataframe.append(df2)
    result = pd.concat(listofdataframe,sort=False)
    print(result.keys())
    print(np.unique(result['interaction_type']))
    print(np.unique(result['label']))
    print("writing to file")
    result.to_hdf('DEFINITIVE.hdf', '/df')

    
    
