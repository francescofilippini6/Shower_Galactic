import numpy as np
import utm
import pandas as pd
import matplotlib.pyplot as plt
import math as m
import sys
import itertools
from scipy.optimize import curve_fit
import seaborn as sns
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
from collections import Counter
from datetime import datetime
from sklearn.metrics import confusion_matrix
import scipy
from matplotlib.colors import LogNorm
import astropy.coordinates as coord
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, concatenate
from astropy.time import Time
import healpy as hp
from skyfield import api
#--------------------------------------------------------
# multi-processor part
#--------------------------------------------------------
import multiprocessing
from functools import partial

def _df_split(tup_arg, **kwargs):
    split_ind, df_split, df_f_name = tup_arg
    return (split_ind, getattr(df_split, df_f_name)(**kwargs))

def df_multi_core(df, df_f_name, subset=None, njobs=-1, **kwargs):
    if njobs == -1:
        njobs = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=njobs)

    try:
        splits = np.array_split(df[subset], njobs)
    except ValueError:
        splits = np.array_split(df, njobs)

    pool_data = [(split_ind, df_split, df_f_name) for split_ind, df_split in enumerate(splits)]
    results = pool.map(partial(_df_split, **kwargs), pool_data)
    pool.close()
    pool.join()
    results = sorted(results, key=lambda x:x[0])
    results = pd.concat([split[1] for split in results])
    return results



#def skyfield_trial(df):
#    ts = api.load.timescale()
#    t = ts.utc(2019, 9, 13, 20)
#    geographic = api.wgs84.latlon(latitude_degrees=42, longitude_degrees=-87)
#    observer = geographic.at(t)
#    pos = observer.from_altaz(alt_degrees=90, az_degrees=0)
#    
#    ra, dec, distance = pos.radec()
    
def reader(filename):
    print("Reading: ",filename)
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
    bdt_cut=0.12
    print("Slicing dataframe at: ", bdt_cut)
    selectedDF=df[df['BDT__cuts_1e2'] > bdt_cut]
    return selectedDF

def weight_astro_spectrum(df):
    spectral_index = 2.0
    spectral_norm = 6 * 10**-12 * 10**(3*spectral_index)  #Neutrino Energy in GeV
    #new_w3_weight=np.array(df['w2'])*spectral_norm*np.power(np.array(np.power(10,df['MCE'])),-spectral_index)
    my_spectral_index = 2.4
    my_spectral_norm = 2.1 * 10**-12 * 10**(3*my_spectral_index)  #Neutrino Energy in GeV
    #imporvements of cpu time of a factor 1000 minimum !!!!
    #new_w3_weight=np.array(df['w2'])*my_spectral_norm*np.power(np.array(np.power(10,df['MCE'])),-my_spectral_index)
    #error_norm=(4.9 *10**-12 * 10**(3*my_spectral_index))*np.power(np.array(np.power(10,df['MCE'])),-my_spectral_index)*np.array(df['w2']) #energy of nu  in GeV
    new_w3_weight=(np.array(df['w2'])*(4.8*10**-7))/(np.power(np.array(np.power(10,df['MCE'])),2.3))*(10**4)*(0.5)
    df['new_w3'] = np.array(new_w3_weight)
    #df['error_norm'] = np.array(error_norm)
    print("done")
    return df


def plot_reco_energy_distro_prediction1(df,dfmuon):
    livetime=3012/365
    fig = plt.figure()
    ax=fig.add_subplot(111)
    #ax.hist(df['TantraEnergy'],bins=50,histtype='step',weights=np.array(df['new_w3'])*livetime,label='nu cosmic total')
    print("first cut")
    df0=df[df['predicted_label']>0.8]
    ax.hist(df0['TantraEnergy'],bins=50,histtype='step',weights=np.array(df0['WeightAtmo'])*livetime,label='nu atmo, prediction 1, evt number:'+ str(sum(df0['WeightAtmo'])))
    df1=df0[df0['label']<0.1]
    ax.hist(df1['TantraEnergy'],bins=50,histtype='step',weights=np.array(df1['WeightAtmo'])*livetime,label='nu atmo, prediction 1, label 0, evt number:'+ str(sum(df1['WeightAtmo'])))
    df1=df0[df0['label']>0.8]
    ax.hist(df1['TantraEnergy'],bins=50,histtype='step',weights=np.array(df1['WeightAtmo'])*livetime,label='nu atmo, prediction 1, label 1, evt number:'+ str(sum(df1['WeightAtmo'])))
    print("second histo")
    ax.hist(df0['TantraEnergy'],bins=50,histtype='step',weights=np.array(df0['new_w3'])*livetime,label='nu cosmic, prediction 1, evt number:'+ str(sum(df0['new_w3'])))
    print("second cut")
    dfmuon=dfmuon[dfmuon['predicted_label']>0.8]
    ax.hist(dfmuon['TantraEnergy'],bins=50,histtype='step',weights=np.array(dfmuon['WeightAtmo']),label='Atmo Muons, prediction 0, evt number:'+ str(sum(dfmuon['WeightAtmo'])))
    ax.set_title('Predicted 1')
    ax.set_xlabel('log10(E/GeV)')
    ax.set_ylabel(r'$\frac{dN}{dE}$')
    ax.set_yscale('log')
    ax.legend()
    plt.show()
    return fig

def confusion_matrix_weighted(df1):
    fig=plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('Total MC sample BDT>0.12')
    print("calculating confusion matrix")
    df_cm=confusion_matrix(df1['label'],df1['predicted_label'],sample_weight=df1['WeightAtmo'])
    df_cm=confusion_matrix(df1['label'],df1['predicted_label'],sample_weight=df1['new_w3'])
    sns.heatmap(df_cm, annot=True,cmap=plt.cm.Blues, fmt='g',ax=ax)
    print("saving figure")
    fig.savefig('Cosmic_weighted_2.3_BDT_0.12_confusion_matrix.png')
    plt.show()
    return fig

#def convergence_angle(lat, lon):
#    """Calculate the converge angle on the UTM grid.
#
#    Parameters
#    ----------
#    lon : number
#        Longitude in rad
#    lat : number
#        Latitude in rad
#
#    """
#    latitude_deg = lat * u.deg
#
#    if latitude_deg > 84 * u.deg or latitude_deg < -80 * u.deg:
#        raise ValueError(
#            "UTM coordinate system is only defined between -80deg S and 84deg N."
#        )
#
#    # detector position, longitude and latitude in rad
#    # lambda  = longitude
#    phi = lat
#
#    # find UTM zone and central meridian
#
#    # longitude of the central meridian of UTM zone in rad
#    lambda0 = longitude_of_central_meridian(utm_zone(lon))
#    omega = lon - lambda0
#    # parameters of the Earth ellipsoid
#    sma = 6378137  # semi-major axis in meters (WGS84)
#    ecc = 0.0066943800  # eccentricity (WGS84)
#
#    rho = sma * (1 - ecc) / pow(1 - ecc * np.sin(phi) ** 2, 3 / 2)
#    nu = sma / np.sqrt(1 - ecc * np.sin(phi) ** 2)
#    psi = nu / rho
#    t = np.tan(phi)
#
#    angle = (
#        np.sin(phi) * omega
#        - np.sin(phi) * omega ** 3 / 3 * pow(np.cos(phi), 2) * (2 * psi ** 2 - psi)
#        - np.sin(phi)
#        * omega ** 5
#        / 15
#        * pow(np.cos(phi), 4)
#        * (
#            psi ** 4 * (11 - 24 * t ** 2)
#            - psi ** 3 * (11 - 36 * t ** 2)
#            + 2 * psi ** 2 * (1 - 7 * t ** 2)
#            + psi * t ** 2
#        )
#        - np.sin(phi)
#        * omega ** 7
#        / 315
#        * pow(np.cos(phi), 6)
#        * (17 - 26 * t ** 2 + 2 * t ** 4)
#    )
#
#    return angle
#
#def utm_zone(lat):
#    """The UTM zone for a given latitude
#
#    Parameters
#    ----------
#    lat : number
#        Latitude in rad
#
#    """
#    return 1 + int((np.pi + lat) / (6 * np.pi / 180))
#
#def longitude_of_central_meridian(utmzone):
#    """The longitude of the central meridian for a given UTM zone.
#    
#    Parameters
#    ----------
#    utmzone : number
#    The UTM zone.
#
#    """
#    zone_width = 6 * np.pi / 180
#    return -np.pi + (utmzone - 1) * zone_width + zone_width / 2
#
def SPEED2coordinateconverter(df):
    antares_northing = 4742381.9
    antares_easting = 268221.6
    antares_height = -2500  # m     (guessed, probably similar to orca)
    antares_utm_zone_number = 32
    antares_utm_zone_letter = "N"
    antares_utm_zone = "{num}{let}".format(num=antares_utm_zone_number, let=antares_utm_zone_letter)
    antares_latitude, antares_longitude = utm.to_latlon(antares_easting, antares_northing, antares_utm_zone_number, antares_utm_zone_letter)
    print(antares_latitude,antares_longitude)
    #antares_lat=42.8  #42?48\' N
    #antares_lon=-6.17  #6?10' E ->  0 / +180 east; 0/-180 west
    locationAntares = locationAntares= EarthLocation(lat=antares_latitude*u.deg , lon= antares_longitude*u.deg, height= antares_height*u.m)
    alt=np.array(np.pi/2-np.arccos(np.array(df['TantraZenith'])))
    azi=np.array(df['TantraAzimuth'])
    print("entering")
    for i in range(10):
        now = datetime.now()
        print(now)
        shift=4.5/24+i*(1-4.5/24-7.2/24)/9
        evt_time=np.array(df['DateMJD'])+shift
        print(evt_time)
        obstime = Time(evt_time,format='mjd')#.to_value('isot')
        frame_hor = AltAz(obstime=obstime, location=locationAntares)
        local_sky=SkyCoord(azi,alt,frame=frame_hor, unit=u.rad)
        print("after local sky")
        #print("Now local to equatorial")
        #equatorial=local_sky.transform_to("icrs")
        print("Now local to galactic")
        gal=local_sky.galactic#transform_to('galactic')
        #print("after transformation")
        #df['eq_ra'] = np.array(equatorial.ra)
        #df['eq_dec'] = np.array(equatorial.dec)
        df['gal_l_offzone'+str(i+1)] = np.array(gal.l)
        df['gal_b_offzone'+str(i+1)] = np.array(gal.b)
        
    print(df.keys())
    print("saving dataframe to file")
    df.to_hdf('Dataframe_zones_astropy.h5', key='df', mode='w')
    return df
    
def angular_uncertanty(df):
    fig = plt.figure()
    df0=df[df['label']<0.1]
    ax=fig.add_subplot(221)
    #ax.plot(df0['MCZenith'],df0['TantraZenith'],'+r')
    ax.hist2d(df0['MCZenith'],df0['TantraZenith'], bins=50)
    ax.set_title('Shower-zenith')
    ax.set_xlabel('mc zenith')
    ax.set_ylabel('Tantra Zenith')
    #Print(df['TantraAngularEstimator'])
    #ax.hist(df['TantraAngularEstimator'],bins=50,histtype='step',label='Tantra Angular estimator distriution')
    #ax.set_title('Tantra Angular Estimator')
    #ax.set_xlabel('Angular Uncertanty')
    #ax.legend()
    df1=df[df['label']<0.8]
    ax2=fig.add_subplot(222)
    ax2.hist2d(df1['MCZenith'],df1['TantraZenith'], bins=50)
    ax2.set_title('Track-zenith')
    ax2.set_xlabel('mc zenith')
    ax2.set_ylabel('Tantra Zenith')
    
    ax3=fig.add_subplot(223)
    #ax.plot(df0['MCZenith'],df0['TantraZenith'],'+r')
    ax3.hist2d(df0['MCAzimuth'],df0['TantraAzimuth'], bins=50)
    ax3.set_title('Shower-Azimuth')
    ax3.set_xlabel('mc azimuth')
    ax3.set_ylabel('Tantra Azimuth')
    
    ax4=fig.add_subplot(224)
    #ax.plot(df0['MCZenith'],df0['TantraZenith'],'+r')
    ax4.hist2d(df1['MCAzimuth'],df1['TantraAzimuth'], bins=50)
    ax4.set_title('Track-Azimuth')
    ax4.set_xlabel('mc azimuth')
    ax4.set_ylabel('Tantra Azimuth')
    
    #ax.plot(df0['MCZenith'],df0['AAZenith'],'+b')
    #ax2.hist(df['TantraAngularEstimator'],bins=50,histtype='step',cumulative=True,label='Cumulative')
    #ax.legend()
    plt.show()
    return fig

def Cumulative_Giulia(df):
    fig = plt.figure()
    df0=df[df['label']<0.1]
    df1=df[df['label']>0.8]
    
    
    ax1=fig.add_subplot(121)
    ax1.hist((180/np.pi)*np.absolute(np.arccos(np.array(df0['TantraZenith']))-np.arccos(np.array(df0['MCZenith']))), bins=1000,label='Tantra',histtype='step', density=True)#,cumulative=True
    ax1.hist((180/np.pi)*np.absolute(np.arccos(np.array(df0['AAZenith']))-np.arccos(np.array(df0['MCZenith']))), bins=1000,label='AA',histtype='step', density=True)#,cumulative=True
    ax1.set_title('not-(a)numuCC')
    ax1.set_xlabel(r'$\Delta\Psi_{zenith}$')
    ax1.set_xlim((0.1,100))
    ax1.set_xscale('log')
    ax1.legend()

    ax2=fig.add_subplot(122)
    ax2.hist((180/np.pi)*np.absolute(np.arccos(np.array(df1['TantraZenith']))-np.arccos(np.array(df1['MCZenith']))), bins=1000,label='Tantra', density=True,histtype='step')#,cumulative=True
    ax2.hist((180/np.pi)*np.absolute(np.arccos(np.array(df1['AAZenith']))-np.arccos(np.array(df1['MCZenith']))), bins=1000,label='AA',density=True,histtype='step')#,cumulative=True
    ax2.set_xlim((0.1,100))
    ax2.set_xscale('log')        
    ax2.set_title('(a)numuCC')
    ax2.set_xlabel(r'$\Delta\Psi_{zenith}$')
    ax2.legend()
    plt.show()
    return fig

def delta_psi_distro(df):
    #replacing nan values with 0 for Mestimator column
    df['Mestimator'] = df['Mestimator'].replace(np.nan, 0)
    print(df['Mestimator'])
    fig = plt.figure()
    df0=df[df['label']<0.1]
    df1=df[df['label']>0.8]
    df00=df0[df0['predicted_label']<0.1]
    df01=df0[df0['predicted_label']>0.8]
    df10=df1[df1['predicted_label']<0.1]
    df11=df1[df1['predicted_label']>0.8]
    
    ax1=fig.add_subplot(221)
    #ax1.hist2d(90-(180/np.pi)*np.arccos(np.array(df00['Mestimator'])),(180/np.pi)*np.absolute(np.arccos(np.array(df00['TantraZenith']))-np.arccos(np.array(df00['MCZenith']))), bins=50,cmap = "RdYlGn_r",weights=df00['WeightAtmo'])
    ax1.hist2d(df00['Mestimator'],(180/np.pi)*np.absolute(np.arccos(np.array(df00['TantraZenith']))-np.arccos(np.array(df00['MCZenith']))), bins=100,cmap = plt.cm.jet,weights=df00['WeightAtmo'],norm=LogNorm())
    ax1.set_title(r'$\Delta\Psi_{zenith}$ Prediction 0, label 0')
    #ax1.set_ylim((0,10))
    ax1.set_xlabel('Mestimator')
    ax1.set_ylabel(r'$\Delta\Psi_{zenith}$')
    
    ax2=fig.add_subplot(222)
    ax2.hist2d(df01['Mestimator'],(180/np.pi)*np.absolute(np.arccos(np.array(df01['TantraZenith']))-np.arccos(np.array(df01['MCZenith']))), bins=100,cmap=plt.cm.jet,norm=LogNorm(),weights=df01['WeightAtmo'])
    ax2.set_title(r'$\Delta\Psi_{zenith}$ Prediction 1, label 0')
    #ax2.set_ylim((0,10))
    #ax2.set_xlabel('mc zenith')
    ax2.set_xlabel('Mestimator')
    ax2.set_ylabel(r'$\Delta\Psi_{zenith}$')
    
    ax3=fig.add_subplot(223)
    ax3.hist2d(df10['Mestimator'],(180/np.pi)*np.absolute(np.arccos(np.array(df10['TantraZenith']))-np.arccos(np.array(df10['MCZenith']))), bins=100,cmap = plt.cm.jet,weights=df10['WeightAtmo'],norm=LogNorm())
    #ax3.plot(90-(180/np.pi)*np.arccos(np.array(df0['MCZenith'])),(180/np.pi)*np.absolute(np.arccos(np.array(df0['TantraZenith']))-np.arccos(np.array(df0['MCZenith']))),'+r')
    ax3.set_title(r'$\Delta\Psi_{zenith}$ Prediction 0, label 1')
    #ax3.set_ylim((0,10))
    #ax3.set_xlabel('mc zenith')
    ax3.set_xlabel('Mestimator')
    ax3.set_ylabel(r'$\Delta\Psi_{zenith}$')
    
    ax4=fig.add_subplot(224) 
    ax4.hist2d(df11['Mestimator'],(180/np.pi)*np.absolute(np.arccos(np.array(df11['TantraZenith']))-np.arccos(np.array(df11['MCZenith']))), bins=100,cmap = plt.cm.jet,weights=df11['WeightAtmo'],norm=LogNorm())
    ax4.set_title(r'$\Delta\Psi_{zenith}$ Predictio 1, label 1')
    #ax4.set_ylim((0,10))
    #ax4.set_xlabel('mc zenith')
    ax4.set_xlabel('Mestimator')
    ax4.set_ylabel(r'$\Delta\Psi_{zenith}$')
    plt.tight_layout()
    
    plt.suptitle("MC zenith and reco uncertainty TANTRA")
    plt.show()
    
    return fig



if __name__ == "__main__":
    filename=sys.argv[1]
    #muon=sys.argv[2]
    df=reader(filename)
    #muon=reader(muon)
    print("LENGTH:",len(df['TriggerT3']))
    #muoncut=cut_dataframe_bdts(muon)
    #df=cut_dataframe_bdts(df)
    print(df.keys())
    #df1=weight_astro_spectrum(df)    
    #Cumulative_Giulia(df)
    #delta_psi_distro(df)
    #angular_uncertanty(df1)
    #original_energy_distro(df1,muoncut)
    #plot_reco_energy_distro_prediction1(df1,muon)
    #plot_reco_energy_distro_prediction0(df1,muon)
    #res = df_multi_core(df=df, df_f_name='isin', njobs=-1, values=lookfor)#subset=['c1'] default None
    SPEED2coordinateconverter(df)
    
    
