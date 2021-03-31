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
from sklearn.metrics import confusion_matrix
import scipy
from matplotlib.colors import LogNorm
import astropy.coordinates as coord
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, concatenate
from astropy.time import Time
import healpy as hp
from skyfield import api
from skyfield.api import wgs84

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
    print("entering")
    evt_time=np.array(df['DateMJD'])
    print(evt_time)
    #azi=np.array(df['TantraAzimuth'])
    azi=np.array(df['MCAzimuth'])
    #true_azimuth = (np.pi / 2 - azi + np.pi + convergence_angle( locationAntares.lat.rad,  locationAntares.lon.rad)) % (2 * np.pi)
    #it is called zenith but it is already the alt!!! Zenith (0,180), alt(-90,90)=90-zenith
    #alt=90-(180/np.pi)*np.arccos(np.array(df['TantraZenith']))
    alt=np.pi/2-np.arccos(np.array(df['MCZenith']))
    obstime = Time(evt_time,format='mjd')#.to_value('isot')
    frame_hor = AltAz(obstime=obstime, location=locationAntares)
    local_sky=SkyCoord(true_azimuth,alt,frame=frame_hor, unit=u.rad)
    print("after local sky")
    print("Now local to equatorial")
    equatorial=local_sky.transform_to("icrs")
    print("Now equatorial to galactic")
    gal=equatorial.transform_to('galactic')
    #print("after transformation")
    df['eq_ra'] = np.array(equatorial.ra)
    df['eq_dec'] = np.array(equatorial.dec)
    df['gal_l'] = np.array(gal.l)
    df['gal_b'] = np.array(gal.b)
    print(df.keys())
    print(len(df['eq_ra']))
    print("saving dataframe to file")
    df.to_hdf('Dataframe_coordinates.h5', key='df', mode='w')
    return df
    
def skyfield_trial(df):
    antares_northing = 4742381.9
    antares_easting = 268221.6
    antares_height = -2500  # m     (guessed, probably similar to orca)
    antares_utm_zone_number = 32
    antares_utm_zone_letter = "N"
    antares_utm_zone = "{num}{let}".format(num=antares_utm_zone_number, let=antares_utm_zone_letter)
    antares_latitude, antares_longitude = utm.to_latlon(antares_easting, antares_northing, antares_utm_zone_number, antares_utm_zone_letter)
    ts = api.load.timescale()
    evt_time=np.array(df['DateMJD'])+ 2400000.5    #from mjd to jd
    #print(evt_time)
    t = ts.tt_jd(evt_time)#now()
    print(len(t))
    alt=(180/np.pi)*(np.pi/2-np.arccos(np.array(df['TantraZenith'])))
    azi=(180/np.pi)*(np.array(df['TantraAzimuth']))
    geographic =wgs84.latlon(latitude_degrees=antares_latitude, longitude_degrees=antares_longitude,elevation_m=antares_height)
    #for a in range(len(evt_time)):
    observer = geographic.at(t[:100])
    print(observer)
    pos = observer.from_altaz(alt_degrees=alt[:100], az_degrees=azi[:100])
    galactic = pos.galactic_latlon()
    print(galactic)
    return 0


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
    #SPEED2coordinateconverter(df)
    skyfield_trial(df)
    
