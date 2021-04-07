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
    #------------
    # dfa ON_9_OFF....
    #-----------
    filename0=sys.argv[1]
    #------------
    # dfa continuos neutrino prediction
    #-----------
    filename1=sys.argv[2]
    #-----------
    dfa=reader(filename0)
    dfb=reader(filename1)
    
    total_label=[]
    for a in np.unique(dfa['interaction_type']):
        dfc=dfb[dfb['interaction_type']==a]
        labels=dfc['predicted_label']
        print(labels)
        dfd=dfa[dfa['interaction_type']==a]
        if len(dfd['label'])==len(labels):
            print(a,len(dfd['label']))
            total_label.append(np.array(labels))
    print(len(total_label))
    total_lab=np.concatenate(total_label, axis=None)
    print(len(total_lab))
    dfa['cont_label_pred']=total_lab
    print("saving to file")
    dfa.to_hdf('ON_OFF_NN_y_continuos.h5', key='df', mode='w')
    
