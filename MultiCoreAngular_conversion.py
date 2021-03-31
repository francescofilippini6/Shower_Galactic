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
#import multiprocessing as mp
#from multiprocesspandas import applyparallel

#import multiprocessing 
#from functools import partial
#def _df_split(tup_arg, **kwargs):
#    split_ind, df_split, df_f_name = tup_arg
#    return (split_ind, getattr(df_split, df_f_name)(**kwargs))
#
#def df_multi_core(df, df_f_name, subset=None, njobs=-1, **kwargs):
#    if njobs == -1:
#        njobs = multiprocessing.cpu_count()
#        pool = multiprocessing.Pool(processes=njobs)
#        
#    try:
#        splits = np.array_split(df[subset], njobs)
#    except ValueError:
#        splits = np.array_split(df, njobs)
#
#    pool_data = [(split_ind, df_split, df_f_name) for split_ind, df_split in enumerate(splits)]
#    results = pool.map(partial(_df_split, **kwargs), pool_data)
#    pool.close()
#    pool.join()
#    results = sorted(results, key=lambda x:x[0])
#    results = pd.concat([split[1] for split in results])
#    return results
#

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
    # if dd  then ,key='df')
    return df

def SPEED2coordinateconverter(dfaa):
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
    Zenith=np.array(dfa['TantraZenith'])
    azi=np.array(dfa['TantraAzimuth'])
    Timemjd=np.array(dfa['DateMJD'])
    print(len(Timemjd))
    alt=np.array(np.pi/2-np.arccos(Zenith))
    #azi=np.array(dfa['TantraAzimuth'])
    print("entering")
    #for i in range(10):
    #now = datetime.now()
    #print(now)
    #shift=4.5/24+i*(1-4.5/24-7.2/24)/9
    evt_time=Timemjd#+shift
    print(evt_time)
    obstime = Time(evt_time,format='mjd')#.to_value('isot')
    frame_hor = AltAz(obstime=obstime, location=locationAntares)
    local_sky=SkyCoord(azi,alt,frame=frame_hor, unit=u.rad)
    print("Now local to galactic")
    gal=local_sky.galactic#transform_to('galactic')
    dfa['gal_l']=gal.l
    dfa['gal_b']=gal.b
    return dfa
    #return dfa.assign(result=(gal.l,gal.b))




if __name__ == "__main__":
    filename=sys.argv[1]
    df=reader(filename)
    print("LENGTH:",len(df['TriggerT3']))
    print(df.keys())
    #res = df_multi_core(df=df, df_f_name='isin', njobs=-1, values=lookfor)#subset=['c1'] default None
    #SPEED2coordinateconverter(df)
    #    df.apply_parallel(SPEED2coordinateconverter, num_processes= (mp.cpu_count() - 1))
    #   df = df.swifter.apply(SPEED2coordinateconverter)
    N=10
    dfp = dd.from_pandas(df, npartitions=N)
    print("Before map partition")
    #print(dfp['DateMJD'].mean().compute())   ok funziona
    #res = dfp.map_partitions(len).compute()   ok funziona
    res = dfp.map_partitions(SPEED2coordinateconverter).compute()
    print("After map partition")
    print(res)

    #result = dfp.compute() 
    
    #print(result)
