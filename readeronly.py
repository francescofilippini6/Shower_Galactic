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

if __name__ == "__main__":
    filename=sys.argv[1]
    df=reader(filename)
    print("LENGTH:",len(df['TriggerT3']))
    if "MUON" in filename:
        print("No label")
    else:
        print(Counter(df['label']))
        print(Counter(df['interaction_type']))
    print(df.keys())
    print(Counter(df['WeightAtmo']))
    #one_labels=[]
    #    for cc in df['TriggerT3']:
    #        if cc==1:
    #            one_labels.append(1)
    #    print(len(one_labels)
    #                  
    # akey=[]
    # for i in df.keys():
    #     for a in df[i]:
    #         if a == -99999:
    #             key.append(i)
    #         else:
    #             continue
    #         
    # print(np.unique(key))
