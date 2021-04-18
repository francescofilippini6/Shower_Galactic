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


def SPEED2coordinateconverter(dfa,cc):
    antares_northing = 4742381.9
    antares_easting = 268221.6
    antares_height = -2500  # m     (guessed, probably similar to orca)
    antares_utm_zone_number = 32
    antares_utm_zone_letter = "N"
    antares_utm_zone = "{num}{let}".format(num=antares_utm_zone_number, let=antares_utm_zone_letter)
    antares_latitude, antares_longitude = utm.to_latlon(antares_easting, antares_northing, antares_utm_zone_number, antares_utm_zone_letter)
    print(antares_latitude,antares_longitude)
    locationAntares = locationAntares= EarthLocation(lat=antares_latitude*u.deg , lon= antares_longitude*u.deg, height= antares_height*u.m)
    Zenith=np.array(dfa.MCZenith)
    print("MCZenith",Zenith)
    print(min(Zenith))
    print(max(Zenith))
    
    azi=np.array(dfa.MCAzimuth)
    alt=np.array(np.pi/2-np.arccos(Zenith))
    Timemjd=np.array(dfa.DateMJD)
    evt_time=Timemjd
    obstime = Time(evt_time,format='mjd')
    frame_hor = AltAz(obstime=obstime, location=locationAntares)
    local_sky=SkyCoord(azi,alt,frame=frame_hor, unit=u.rad)
    print("Now local to galactic")
    gal=local_sky.galactic
    return (gal.l.deg,gal.b.deg)

def trigger_counter_plot(df):
    fig=plt.figure()
    ax=fig.add_subplot(121)
    ax.hist(df['TriggerCounter'],bins=50,histtype='step')
    ax.set_title("Trigger counter distro")
    ax.set_xlabel('TriggerCounter value')
    ax.set_ylabel(r'entries')
    ax.set_yscale('log')
    ax1=fig.add_subplot(122)
    ax1.plot(np.power(10,df['MCE']),df['TriggerCounter'],'+r')
    ax1.set_title("Trigger counter vs MCE")
    ax1.set_xlabel('MCE')
    ax1.set_ylabel('Trigger Counter')
    ax1.set_yscale('log')
    plt.show()
    return fig

if __name__ == "__main__":
    filename=sys.argv[1]
    df=reader(filename)
    
    print("LENGTH:",len(df['TantraEnergy']))
    #    if "MUON" in filename:
    #        print("No label")
    #    else:
    #        print(Counter(df['label']))
    #        print(Counter(df['interaction_type']))

    print(df.keys())
    
    print(Counter(df['label']))
    #print(Counter(df['interaction_type']))
    print('3N:',Counter(df['Trigger3N']))
    print('T3:',Counter(df['TriggerT3']))
    #print(Counter(df['TriggerCounter']))

    #trigger_counter_plot(df)

    # #df=df.head(10)
    # print("usual conversion",SPEED2coordinateconverter(df,1))
    # ra=np.array(df['MCRa'])
    # dec=np.array(df['MCDec'])
    # print("RA",ra)
    # print("DEC",dec)
    # asky=SkyCoord(ra=ra*u.rad,dec=dec*u.rad,frame='icrs')
    # ga=asky.galactic
    # print(ga.l.deg)
    # print(ga.b.deg)
    # #print(df['predicted_label'])
    #    df1=df[df['interaction_type']=='numuCC']
    #    df2=df[df['interaction_type']=='numuNC']
    #
    #    df1= df1[['Mestimator', 'TantraZenith', 'TantraAzimuth', 'Lambda', 'Beta',
    #       'WeightAtmo', 'TantraEnergy', 'RunDurationYear', 'DateMJD', 'MCE',
    #       'MCZenith', 'MCAzimuth', 'MCX', 'MCY', 'MCZ', 'MCRho', 'w2', 'w3',
    #       'MCRa', 'MCDec', 'BDT_default', 'BDT__cuts_1e2', 'TantraRho']]
    #
    #    df2= df2[['Mestimator', 'TantraZenith', 'TantraAzimuth', 'Lambda', 'Beta',
    #       'WeightAtmo', 'TantraEnergy', 'RunDurationYear', 'DateMJD', 'MCE',
    #       'MCZenith', 'MCAzimuth', 'MCX', 'MCY', 'MCZ', 'MCRho', 'w2', 'w3',
    #       'MCRa', 'MCDec', 'BDT_default', 'BDT__cuts_1e2', 'TantraRho']]
    #    
    #print(df1.head(10))
    #print(df2.head(10))
    
    #print(Counter(df['WeightAtmo']))
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
