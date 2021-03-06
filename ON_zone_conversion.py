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
#import healpy as hp
#import swifter
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
    print("CC:",cc)
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
    Zenith=np.array(dfa.TantraZenith)
    azi=np.array(dfa.TantraAzimuth)
    print("azi",azi)
    Timemjd=np.array(dfa.DateMJD)
    print("Array:",Timemjd)
    print("LEN:",len(Timemjd))
    print("TYPE:",Timemjd.dtype)
    alt=np.array(np.pi/2-np.arccos(Zenith))
    #shift=4.5/24+i*(1-4.5/24-7.2/24)/9
    evt_time=Timemjd#+shift
    #print(evt_time)
    obstime = Time(evt_time,format='mjd')#.to_value('isot')
    frame_hor = AltAz(obstime=obstime, location=locationAntares)
    local_sky=SkyCoord(azi,alt,frame=frame_hor, unit=u.rad)
    print("Now local to galactic")
    gal=local_sky.galactic#transform_to('galactic')
    dfa['gal_l']=gal.l.deg
    dfa['gal_b']=gal.b.deg
    #dfa.assign(gal_lat=gal.l,gal_lon=gal.b)
    #return dfa
    #
    return dfa
     



if __name__ == "__main__":
    filename=sys.argv[1]
    df=reader(filename)
    print("LENGTH:",len(df['TriggerT3']))
    print(df.keys())
    N=8
    dfp = dd.from_pandas(df, npartitions=N)
    print("Before map partition")
    res = dfp.map_partitions(SPEED2coordinateconverter,1)
    dfi=res.compute()
    print("After map partition")
    #print(dfi)
    #print(dfi['TantraZenith'],dfi['TantraAzimuth'],dfi['DateMJD'],dfi['gal_b'],dfi['gal_l'])
    print("writing file")
    #dfi.to_hdf('Data_converted_partition.h5', key='df', mode='w')
    dfi.to_hdf('ONzone_converted.hdf', '/df')
