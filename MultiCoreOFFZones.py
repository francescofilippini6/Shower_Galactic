import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math as m
import sys
import time
import itertools
import astropy.coordinates as coord
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, concatenate
from astropy.time import Time
import dask.dataframe as dd
import utm
from dask.multiprocessing import get
#--------------------------------------------------------
# multi-processor part
#--------------------------------------------------------
def reader(filename):
    print("Reading: ",filename)
    df = pd.read_hdf(filename)
    # if dd  then ,key='df')
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
    Zenith=np.array(dfa.TantraZenith)
    azi=np.array(dfa.TantraAzimuth)
    alt=np.array(np.pi/2-np.arccos(Zenith))
    Timemjd=np.array(dfa.DateMJD)
    for i in range(10):   #for i in range(10):
        shift=4.5/24+i*(1-4.5/24-7.2/24)/9
        evt_time=Timemjd+shift
        obstime = Time(evt_time,format='mjd')
        frame_hor = AltAz(obstime=obstime, location=locationAntares)
        #look at the order of insertion azi e alt
        local_sky=SkyCoord(az=azi,alt=alt,frame=frame_hor, unit=u.rad)
        print("Now local to galactic")
        gal=local_sky.galactic
        dfa['gal_l'+str(i)]=gal.l.deg
        dfa['gal_b'+str(i)]=gal.b.deg
        
    return dfa
     



if __name__ == "__main__":
    filename=sys.argv[1]
    df=reader(filename)
    #print("LENGTH:",len(df['TriggerT3']))
    df=df.drop(['AAZenith', 'AAAzimuth','Trigger3N', 'TriggerT3','IntegralCharge','MeanCharge', 'StdCharge','TriggerCounter','TantraLines','TantraHits', 'Trigger3N', 'TriggerT3', 'Nrun','TriggerCounter', 'FrameIndex', 'MCNGenEvent', 'MCNuType','MCLeadingParticle', 'MCMultiplicity',  'TantraAngularEstimator', 'TantraX', 'TantraY','TantraZ', 'NOnTime', 'GridQuality', 'TrackLength','WeightAstro','NEW_LIKELIHOOD_3D_ATMO', 'NEW_LIKELIHOOD_3D_ASTRO'],axis=1)
    print(df.keys())
    N=8
    #df=df.head(100)
    dfp = dd.from_pandas(df, npartitions=N)
    print("Before map partition")
    start = time.process_time()
    res = dfp.map_partitions(SPEED2coordinateconverter,1)
    dfi=res.compute()
    print("EXECUTION TIME:",time.process_time() - start)
    
    print("writing file")
    #dfi.to_hdf('ON_OFF_converted_partition.h5', key='df', mode='w')
    dfi.to_hdf('OFFZones_Muon_bdtless0.12.hdf', '/df')    
