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
from scipy import special
from scipy.stats import poisson
import scipy.stats as st
import scipy.special as sc

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

def weight_astro_spectrum(df):
    spectral_index = 2.0
    spectral_norm = 6 * 10**-12 * 10**(3*spectral_index)  #Neutrino Energy in GeV
    #new_w3_weight=np.array(df['w2'])*spectral_norm*np.power(np.array(np.power(10,df['MCE'])),-spectral_index)
    my_spectral_index = 2.5
    my_spectral_norm = 3.5 * 10**-12 * 10**(3*my_spectral_index)  #Neutrino Energy in GeV
    #imporvements of cpu time of a factor 1000 minimum !!!!
    #new_w3_weight=np.array(df['w2'])*my_spectral_norm*np.power(np.array(np.power(10,df['MCE'])),-my_spectral_index)
    #error_norm=(4.9 *10**-12 * 10**(3*my_spectral_index))*np.power(np.array(np.power(10,df['MCE'])),-my_spectral_index)*np.array(df['w2']) #energy of nu  in GeV
    new_w3_weight=(np.array(df['w2'])*(4.8*10**-7))/(np.power(np.array(np.power(10,df['MCE'])),2.3))*(10**4)*(0.5)
    df['new_w3'] = np.array(new_w3_weight)
    #df['error_norm'] = np.array(error_norm)
    print("done")
    return df


def new_weighter(df):
    print("inside new weighter")
    identifier=list(range(len(df['interaction_type'])))
    df['identifier']=identifier
    weight=[]
    #excluding muons
    print("first slice, muon out")
    dfslice=df[df['interaction_type'].str.contains("nu")]
    #excluding particles out of the direction
    dfslice=zonesSelection(dfslice,'MC_gal_b', 'MC_gal_l')
    new_w3_weight=(np.array(dfslice['w2'])*(4.8*10**-7))/(np.power(np.array(np.power(10,dfslice['MCE'])),2.3))*(10**4)*(0.5)/0.023
    dfslice['new_w3'] = np.array(new_w3_weight)
    print("done")
    #df=df.reset_index().melt(id_vars='gene', value_name='quantile')
    print("merging")
    final=pd.merge(df,dfslice[['identifier','new_w3']], on='identifier', how='left')
    #final=df.merge(dfslice, on='identifier', how='left')
    final['new_w3'] = final['new_w3'].fillna(0)
    final=final.drop(['identifier'],axis=1)
    print(final.keys())
    print("out")
    return final

def zonesSelection(df,lat,lon):
    print("cut on b latitude")
    df=df[np.absolute(df[lat])<3]
    print("cut on l longitude")
    df=df[(df[lon]<40) | (df[lon]>320)]
    return df

def converterMC(df):
    antares_northing = 4742381.9
    antares_easting = 268221.6
    antares_height = -2500  # m     (guessed, probably similar to orca)
    antares_utm_zone_number = 32
    antares_utm_zone_letter = "N"
    antares_utm_zone = "{num}{let}".format(num=antares_utm_zone_number, let=antares_utm_zone_letter)
    antares_latitude, antares_longitude = utm.to_latlon(antares_easting, antares_northing, antares_utm_zone_number, antares_utm_zone_letter)
    print(antares_latitude,antares_longitude)
    locationAntares = locationAntares= EarthLocation(lat=antares_latitude*u.deg , lon= antares_longitude*u.deg, height= antares_height*u.m)
    Zenith=np.array(df['MCZenith'])
    azi=np.array(df['MCAzimuth'])
    alt=np.array(np.pi/2-np.arccos(Zenith))
    Timemjd=np.array(df['DateMJD'])
    obstime = Time(Timemjd,format='mjd')
    frame_hor = AltAz(obstime=obstime, location=locationAntares)
    local_sky=SkyCoord(azi,alt,frame=frame_hor, unit=u.rad)
    print("Now local to galactic")
    gal=local_sky.galactic
    df['gal_l_MC']=gal.l.deg
    df['gal_b_MC']=gal.b.deg
    return df

def plot_ONZONE(df,bins):
    #livetime=3012/365
    #fig = plt.figure()
    #ax=fig.add_subplot(111)
    print("ON Zone")
    #---------------------------------------------
    # cut on Tantra directions to be in ON
    #---------------------------------------------
    
    df=zonesSelection(df,'gal_b','gal_l')
    print("ON histo")
    hist0, _ = np.histogram(df['TantraEnergy'], bins=bins,weights=np.array(df['WeightAtmo']))
    #print("OFF events:", sum(hist0))
    #--------------------------------------------------
    #  also MC zenith and azimuth must be in the ON region
    #--------------------------------------------------
    #df=converterMC(df)
    #df=zonesSelection(df,'gal_b_MC','gal_l_MC')

    hist1, _ = np.histogram(df['TantraEnergy'], bins=bins,weights=np.array(df['new_w3']))    
    
    print("Cosmic signal events:", sum(hist1))
    print("ON events atm:", sum(hist0))
    histsumON=hist0+hist1
    print("ON events tot:", sum(histsumON))
    return (histsumON,hist1)


def plot_OFFZONE(df):
    print("OFF_Zones_selection")
    df0=zonesSelection(df,'gal_b0','gal_l0')
    df1=zonesSelection(df,'gal_b1','gal_l1')
    df2=zonesSelection(df,'gal_b2','gal_l2')
    df3=zonesSelection(df,'gal_b3','gal_l3')
    df4=zonesSelection(df,'gal_b4','gal_l4')
    df5=zonesSelection(df,'gal_b5','gal_l5')
    df6=zonesSelection(df,'gal_b6','gal_l6')
    df7=zonesSelection(df,'gal_b7','gal_l7')
    df8=zonesSelection(df,'gal_b8','gal_l8')
    print("OFF histo")
    
    hist0, bins = np.histogram(df0['TantraEnergy'], bins=50,weights=np.array(df0['WeightAtmo']))
    hist1, _ = np.histogram(df1['TantraEnergy'], bins=bins,weights=np.array(df1['WeightAtmo']))
    hist2, _ = np.histogram(df2['TantraEnergy'], bins=bins,weights=np.array(df2['WeightAtmo']))
    hist3, _ = np.histogram(df3['TantraEnergy'], bins=bins,weights=np.array(df3['WeightAtmo']))
    hist4, _ = np.histogram(df4['TantraEnergy'], bins=bins,weights=np.array(df4['WeightAtmo']))
    hist5, _ = np.histogram(df5['TantraEnergy'], bins=bins,weights=np.array(df5['WeightAtmo']))
    hist6, _ = np.histogram(df6['TantraEnergy'], bins=bins,weights=np.array(df6['WeightAtmo']))
    hist7, _ = np.histogram(df7['TantraEnergy'], bins=bins,weights=np.array(df7['WeightAtmo']))
    hist8, _ = np.histogram(df8['TantraEnergy'], bins=bins,weights=np.array(df8['WeightAtmo']))
    
    histsum=(hist0+hist1+hist2+hist3+hist4+hist5+hist6+hist7+hist8)/9
    print("mean OFF events:",sum(histsum))
    print(bins)
    return (histsum,bins)

def ONandOFF(df):
    offhisto,binning=plot_OFFZONE(df)
    onhisto,cosmic=plot_ONZONE(df,binning)
    center = (binning[:-1] + binning[1:]) / 2
    cumulativeon=[]
    cumulativeoff=[]
    cumulativecosmic=[]
    for i in range(len(binning)-1):
        cumulativeon.append(sum(onhisto[i:]))
        cumulativeoff.append(sum(offhisto[i:]))
        cumulativecosmic.append(sum(cosmic[i:]))

    saver = pd.DataFrame({'OnEntries': onhisto,
                          'OffEntries': offhisto,
                          'cosmic' : cosmic,
                          'center_bin' : center,
                          'bin_edges': binning[1:]})
    
    #saver.to_csv('Energy_Mrf.csv')
    #------------------------------------------------------
    # relative discrepancy 
    #------------------------------------------------------
    #discrepancy=(np.absolute(onhisto-offhisto)/offhisto)*100

    #------------------------------------------------------
    # ratio discrepancy 
    #------------------------------------------------------
    #discrepancy=onhisto/offhisto
    #------------------------------------------------------
    # sigma discrepancy 
    #------------------------------------------------------

    discrepancy=[]
    for binc in range(len(onhisto)):
        tailpos=poisson.cdf(onhisto[binc],offhisto[binc]) #cdf=cumulative distribution function
        pvalue=1-tailpos
        discrepancy.append(st.norm.ppf(1-pvalue))

    #discrepancy1=[]
   # for bind in range(len(onhisto)):
   #     LieMa=np.sqrt(2)*np.sqrt(onhisto[bind]*np.log(2*onhisto[bind]/(onhisto[bind]+offhisto[bind]))+offhisto[bind]*np.log(2*offhisto[bind]/(onhisto[bind]+offhisto[bind])))
   #     discrepancy1.append(LieMa)
        
    #discrepancy2=[]
    #for bine in range(len(onhisto)):
    #    onoffarticle=sc.betainc(0.5,onhisto[bine],1+offhisto[bine])/sc.beta(onhisto[bine],1+offhisto[bine])
    #    discrepancy2.append(onoffarticle)
    #print(discrepancy2)
    
    #print(discrepancy,len(discrepancy))    
    fig = plt.figure()
    gs = fig.add_gridspec(nrows=4, ncols=2,  hspace=0)#gs=GridSpec(4,2)
    ax=fig.add_subplot(gs[:-1,0])
    ax.plot(center,offhisto,'+r',label='offzone')
    ax.plot(center,onhisto,'+b',label='onzone')
    ax.plot(center,cosmic,'+g',label='onzone cosmic')
    ax.axvline(x=3.582,color='k',linestyle='--')
    ax.set_xlabel('log10(E/GeV)')
    ax.set_ylabel(r'$\frac{dN}{dE}$')
    ax.set_yscale('log')
    ax.legend()
    ax1=fig.add_subplot(122)
    ax1=fig.add_subplot(gs[:,1])
    ax1.plot(center,cumulativeoff,'--r',label='offzone')
    ax1.plot(center,cumulativeon,'--b',label='onzone')
    ax1.plot(center,cumulativecosmic,'--g',label='onzone cosmic')
    
    ax1.set_xlabel('log10(E/GeV)')
    ax1.set_title('Cumulative')
    ax1.set_ylabel(r'$\frac{dN}{dE}$')
    ax1.set_yscale('log')
    ax1.legend()
    #x=np.arange(min(binning),max(binning),1000)
    ax2=fig.add_subplot(gs[-1:,0])
    #ax2.fill_between(x, 0, 2,color='green', alpha=0.5)
    #ax2.fill_between(x, 2, 3,color='yellow', alpha=0.5)
    #ax2.fill_between(x, 3, 10,color='red', alpha=0.5)
    ax2.axhline(y=2, color='g', linestyle='-',label=r'2 $\sigma$')
    ax2.axhline(y=3, color='r', linestyle='-',label=r'3 $\sigma$')
    ax2.axhline(y=5, color='b', linestyle='-',label=r'5 $\sigma$')
    ax2.plot(center,discrepancy,'+k',label=r'poisson fluctuation ($\mu_s$=0)')
    #ax2.plot(center,discrepancy1,'xb',label='Li and Ma')
    #ax2.plot(center,discrepancy2,'.r',label='Cousins Linnemann Tucker')
    ax2.legend()
    
    ax2.set_ylabel(r'$\sigma$')
    ax2.set_xlabel('log10(E/GeV)')
    ax2.set_yscale('log')
    plt.show()
    return (offhisto,onhisto,cosmic)
    

if __name__ == "__main__":
    listofdataframe=[]
    #------------------------------------------------
    filename=sys.argv[1]
    df=reader(filename)
    #df=df.head(10)
    #df1=df.drop(['predicted_label'],axis=1)
    #labels=np.array(df1['cont_label_pred'])
    #df1=df1.drop(['cont_label_pred'],axis=1)
    #df1['predicted_label']=labels
    df=df.drop(['Mestimator', 'TantraZenith', 'TantraAzimuth', 'Lambda', 'Beta', 'RunDurationYear', 'MCX', 'MCY', 'MCZ', 'MCRho', 'w3','MCRa', 'MCDec', 'BDT_default', 'TantraRho', 'label'],axis=1)
    
    #print(df.keys())
    muon=sys.argv[2]
    dfm=reader(muon)
    #dfm=dfm.head(10)
    dfm=dfm.drop(['Mestimator', 'TantraZenith', 'TantraAzimuth', 'Lambda', 'Beta', 'RunDurationYear', 'MCX', 'MCY', 'MCZ', 'MCRho', 'w3','MCRa', 'MCDec', 'BDT_default', 'TantraRho', 'label'],axis=1)
    
    #df1=weight_astro_spectrum(df)
    df=new_weighter(df)
    dfm=new_weighter(dfm)
    listofdataframe.append(df)
    listofdataframe.append(dfm)
    dff = pd.concat(listofdataframe,sort=False)
    dfbdt=dff[dff['BDT__cuts_1e2']>0.15]
    dfnn=dfbdt[dfbdt['predicted_dropout']<0.3]
    print(dff.keys())
    OFF,ON,cosmic=ONandOFF(dfnn)
    


    
