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
    df = pd.read_hdf(filename)
    return df


def cut_dataframe_bdts(df):
    print("Bdt cut")
    bdt_cut=0.10
    selectedDF=df[df['BDT__cuts_1e2'] > bdt_cut]
    return selectedDF
 
def tantrarhoinsertion(df):
    x=np.array(df['TantraX'])
    y=np.array(df['TantraY'])
    rho=np.sqrt(x**2+y**2)
    return rho

if __name__ == "__main__":
    filename=sys.argv[1].split(',')
    listofdf=[]
    for a in filename:
        print('retrieving file: ',a)
        df=reader(a)
        #df1=cut_dataframe_bdts(df)
        listofdf.append(df)
    print(len(listofdf))
    result = pd.concat(listofdf)
    rh=tantrarhoinsertion(result)
    print(len(result['TantraX']))
    print(len(rh))
    result['TantraRho']=rh
    result['TantraEnergy']=np.log10(np.array(result['TantraEnergy']))
    result['TantraZenith']=np.cos(np.array(result['TantraZenith']))
    result['AAZenith']=np.cos(np.array(result['AAZenith']))
    #print('results',len(result['BDT__cuts_1e2']))
    print('results',result.keys)
    #print(result.keys)
    result.to_hdf('ShowerData_BDT_greater_0.10.h5', key='df', mode='w')
    
