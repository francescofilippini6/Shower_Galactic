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
    bdt_cut=0.33
    selectedDF=df[df['BDT__cuts_1e2'] > bdt_cut]
    return selectedDF
 


if __name__ == "__main__":
    filename=sys.argv[1].split(',')
    listofdf=[]
    for a in filename:
        print('retrieving file: ',a)
        df=reader(a)
        df1=cut_dataframe_bdts(df)
        listofdf.append(df1)
    print(len(listofdf))
    result = pd.concat(listofdf)
    print(result.keys)
    print('results',len(result['BDT__cuts_1e2']))
    result.to_hdf('ShowerData_BDT_greater_0.10.h5', key='df', mode='w')
    
