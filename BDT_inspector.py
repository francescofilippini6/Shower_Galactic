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


if __name__ == "__main__":
    filename=sys.argv[1]
    df=reader(filename)
    
    print("LENGTH:",len(df['TantraEnergy']))
    
    print(df.keys())
    
    print(Counter(df['label']))
    print(Counter(df['interaction_type']))
    aa=np.array(df['BDT__cuts_1e2'])
    print("BDT min",min(aa))
    print("BDT max",max(aa))
    
