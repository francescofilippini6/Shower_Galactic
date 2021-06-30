import numpy as np
from matplotlib import gridspec
import  csv
import pandas as pd
import matplotlib.pyplot as plt
import math as m
import utm
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


def plotter(dfdata,dfmc,key):
    print(key)
    #bins=0
    datahist, bins = np.histogram(dfdata[key],bins=40)
    mchist, _ = np.histogram(dfmc[key],bins=bins,weights=np.array(dfmc['WeightAtmo'])+np.array(dfmc['']))
    print(bins)
    print(mchist)
    #center = (bins[:-1] + bins[1:]) / 2
    discrepancy=(datahist/mchist)
    #discrepancy=[]
    #for u in mchist:
    #    if u==0:
    #        discrepancy.append(0)
    #    else:
    #        discrepancy.append(np.absolute(mchist-datahist)/sqrt(mchist))
    
    fig = plt.figure()
    gs = fig.add_gridspec(nrows=4, ncols=1,  hspace=0)#gs=GridSpec(4,2)
    ax=fig.add_subplot(gs[:-1,0])
    #ax=fig.add_subplot(111)
    ax.set_title(key)
    hh , hbin , aaa = ax.hist(bins[:-1], bins, weights=mchist,histtype='step',label='MC')
    ax.hist(bins[:-1], bins, weights=datahist,histtype='step',label='DATA')
    ax.legend()
    ax.set_yscale('log')
    center = (hbin[:-1] + hbin[1:]) / 2
    # ax.plot(center,mchist,'+r',label='MC')
    #ax.plot(center,datahist,'+b',label='DATA')
    print(center)
    ax2=fig.add_subplot(gs[-1:,0])
    ax2.plot(center,discrepancy,'+k')
    ax2.set_ylabel('data/mc')
    plt.savefig('/Users/francescofilippini/Desktop/Shower_event_analysis/MC_dat_comparison/'+key+'.png')
    return fig

if __name__ == "__main__":
    #------------------------------------------------
    MC=sys.argv[1]
    dfMonC=reader(MC)
    #dfMonC=dfMonC[dfMonC['BDT__cuts_1e2']>0.10]
    print('1')
    data=sys.argv[2]
    dfdata=reader(data)
    print('2')
    muon=sys.argv[3]
    dfmuon=reader(muon)
    dfmuon=dfmuon[dfmuon['BDT__cuts_1e2']>0.10]
    print('2.5')
    listdf=[]
    listdf.append(dfMonC)
    listdf.append(dfmuon)
    dfmc = pd.concat(listdf,sort=False)
    print('2.5')
    #dfmc=dfmc[dfmc['BDT__cuts_1e2']>0.10]
    dfdata=dfdata[dfdata['BDT__cuts_1e2']>0.10]
    print('3')
    listofkey=['TantraLines', 'TantraHits', 'Mestimator', 'TantraZenith','TantraAzimuth', 'TantraAngularEstimator', 'TantraX', 'TantraY','TantraZ', 'Lambda','Beta', 'TrackLength','TantraEnergy','TantraRho','TriggerCounter','NOnTime','AAZenith', 'AAAzimuth','GridQuality','Trigger3N', 'TriggerT3','NEW_LIKELIHOOD_3D_ATMO','IntegralCharge','MeanCharge', 'StdCharge']
    key1=['predicted_dropout']
    #key1=['predicted_label']
    for keyy in key1:
        plotter(dfdata,dfmc,keyy)
        
        
