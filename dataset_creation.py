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
    selectedDF=df[df['BDT__cuts_1e2'] > bdt_cut]
    return selectedDF

def correlation_matrix(df):
    #df=df.head(100)
    df1 = df[['TantraLines', 'TantraHits', 'Mestimator', 'TantraZenith',
              'TantraAzimuth', 'TantraAngularEstimator', 'TantraX', 'TantraY',
              'TantraZ', 'Lambda','Beta', 'TrackLength','TantraEnergy','TantraRho','TriggerCounter','NOnTime','AAZenith', 'AAAzimuth','GridQuality','Trigger3N', 'TriggerT3']]
    corrMatrix = df1.corr()
    sn.heatmap(corrMatrix)
    plt.show()

def column_selector(df,name):
    #df1 = df[['TantraLines', 'TantraHits', 'Mestimator', 'TantraZenith',
    #          'TantraAzimuth', 'TantraAngularEstimator', 'TantraX', 'TantraY',
    #          'TantraZ','Lambda','Beta', 'TrackLength','TantraEnergy','TantraRho','IntegralCharge','MeanCharge', 'StdCharge',
    #          'TriggerCounter','GridQuality','AAZenith', 'AAAzimuth','Trigger3N', 'TriggerT3','NOnTime','MCE','MCZenith','MCAzimuth']]
    y=[]
    if 'numuCC' in name:
        y=np.ones(len(df1['TantraLines']))

    else:
        y=np.zeros(len(df1['TantraLines']))
        #y=[0]*len(df1['TantraLines'])
    df1['label']=y
    return df1

def type_appender(df,name):
    interaction_type=name.split('.')[0]
    print(interaction_type)
    y=[]
    for a in range(len(df['TantraLines'])):
        y.append(interaction_type)
    df['interaction_type']=y
    return df

def plotter(df1,df2,name):
    print("inside function")
    interesting=['NOnTime','TantraHits','TriggerCounter','TantraEnergy','Mestimator','IntegralCharge']
    fig, ((ax1, ax2,ax5), (ax3, ax4,ax6)) = plt.subplots(2, 3)
    fig.suptitle(name, fontsize=12)
    ax1.set_title('NOnTime')
    ax1.set_yscale('log')
    df1['NOnTime'].hist(bins=100,alpha=0.5,ax=ax1,label='BDT_cut')#,density=True)
    df2['NOnTime'].hist(bins=100,alpha=0.5,ax=ax1,label='random_sample')#,density=True)
    ax1.legend()
    ax2.set_title('TantraHits')
    ax2.set_yscale('log')
    df1['TantraHits'].hist(bins=100,alpha=0.5,ax=ax2,label='BDT_cut')#,density=True)
    df2['TantraHits'].hist(bins=100,alpha=0.5,ax=ax2,label='random_sample')#,density=True)
    ax2.legend()
    ax3.set_title('TriggerCounter')
    ax3.set_yscale('log')
    df1['TriggerCounter'].hist(bins=100,alpha=0.5,ax=ax3,label='BDT_cut')#,density=True)
    df2['TriggerCounter'].hist(bins=100,alpha=0.5,ax=ax3,label='random_sample')#,density=True)
    ax3.legend()
    ax4.set_title('TantraEnergy')
    ax4.set_yscale('log')
    df1['TantraEnergy'].hist(bins=100,alpha=0.5,ax=ax4,label='BDT_cut')#,density=True)
    df2['TantraEnergy'].hist(bins=100,alpha=0.5,ax=ax4,label='random_sample')#,density=True)
    ax4.legend()
    ax5.set_title('Mestimator')
    ax5.set_yscale('log')
    df1['Mestimator'].hist(bins=100,alpha=0.5,ax=ax5,label='BDT_cut')#,density=True)
    df2['Mestimator'].hist(bins=100,alpha=0.5,ax=ax5,label='random_sample')#,density=True)
    ax5.legend()
    ax6.set_title('IntegralCharge')
    ax6.set_yscale('log')
    df1['IntegralCharge'].hist(bins=100,alpha=0.5,ax=ax6,label='BDT_cut')#,density=True)
    df2['IntegralCharge'].hist(bins=100,alpha=0.5,ax=ax6,label='random_sample')#,density=True)
    ax6.legend()
    plt.show()

def aa_3D_histogram(listofdf):
    fig=plt.figure()
    df=pd.concat(listofdf[:6])
    df1=pd.concat(listofdf[-2:])
    ax = fig.add_subplot(111)
    ax.scatter(df['NOnTime'],df['TantraHits'])
    ax.scatter(df1['NOnTime'],df1['TantraHits'])
    plt.show()

#def directions(df):
#    fig=plt.figure()
#    ax = fig.add_subplot(231)
#    ax.set_title('Azimuth')
#    ax.set_xlabel('TANTRA')
#    ax.set_xlabel('AA')
#    ax.scatter(df['TantraAzimuth'],df['AAAzimuth'])
#    ax1 = fig.add_subplot(232)
#    ax1.scatter(df1['MCAzimuth'],df1['AAAzimuth'])
#    ax1.set_title('Azimuth')
#    ax1.set_xlabel('MC')
#    ax1.set_xlabel('AA')
#    ax1 = fig.add_subplot(233)
#    ax1.scatter(df1['MCAzimuth'],df1['TantraAzimuth'])
#    ax1.set_title('Azimuth')
#    ax1.set_xlabel('TANTRA')
#    ax1.set_xlabel('AA')
#    ax1 = fig.add_subplot(232)
#    ax1.scatter(df1['TantraZenith'],df1['AAAZenith'])
#    ax1.set_title('Zenith')
#    ax1.set_xlabel('TANTRA')
#    ax1.set_xlabel('AA')
#    ax1 = fig.add_subplot(232)
#    ax1.scatter(df1['TantraZenith'],df1['AAAZenith'])
#    ax1.set_title('Zenith')
#    ax1.set_xlabel('TANTRA')
#    ax1.set_xlabel('AA')
#    ax1 = fig.add_subplot(232)
#    ax1.scatter(df1['TantraZenith'],df1['AAAZenith'])
#    ax1.set_title('Zenith')
#    ax1.set_xlabel('TANTRA')
#    ax1.set_xlabel('AA')
   
if __name__ == "__main__":
    filename=sys.argv[1].split(',')
    listofdataframe=[]
    #datasetlength=400000
    for counter,a in enumerate(filename):
        print(a)
        df=reader(a)
        print("Before cut",len(df['TantraX']))
        df1=cut_dataframe_bdts(df)
        print("after cut",len(df1['TantraX']))
        df2=column_selector(df1,a)
        df3=type_appender(df2,a)
        listofdataframe.append(df3)
        #if 'numuCC' in a:
        #listofdataframe.append(column_selector(df,a))
        #df2=df1.sample(n = 1130000)
        #aa_3D_histogram(df2)
        #listofdataframe.append(column_selector(df2,a))
        #plotter(df1,df2,a)
        #else:
            #listofdataframe.append(column_selector(df,a))
            #df2=df1.sample(n = 380000)
            #aa_3D_histogram(df2)
            #listofdataframe.append(column_selector(df2,a))
            #plotter(df1,df2,a)
            #plotter(listofdataframe[counter])
    #aa_3D_histogram(listofdataframe)
    result = pd.concat(listofdataframe)
    print(result.keys())
    print(np.unique(result['interaction_type']))
    print(np.unique(result['label']))
    #print(result)
    #print("Final dataset length",len(result['label']))
    result.to_hdf('MC_all_dataset_BDT_CUT_0.12.h5', key='df', mode='w')
    #coordinate_plooter(df)
    #correlation_matrix(df)

