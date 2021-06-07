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

def SPEED2coordinateconverter():
    antares_northing = 4742381.9
    antares_easting = 268221.6
    antares_height = -2500  # m     (guessed, probably similar to orca)
    antares_utm_zone_number = 32
    antares_utm_zone_letter = "N"
    antares_utm_zone = "{num}{let}".format(num=antares_utm_zone_number, let=antares_utm_zone_letter)
    antares_latitude, antares_longitude = utm.to_latlon(antares_easting, antares_northing, antares_utm_zone_number, antares_utm_zone_letter)
    print(antares_latitude,antares_longitude)
    locationAntares = locationAntares= EarthLocation(lat=antares_latitude*u.deg , lon= antares_longitude*u.deg, height= antares_height*u.m)
    zenith=[-0.628058,-0.895332,-0.118879]
    azimuth=[4.937221,3.810610,0.785857]
    timem=[56543.602280,56543.674005,56543.594699]
    azi=np.array(azimuth)
    alt=np.array(np.pi/2-np.arccos(np.array(zenith)))
    evt_time=np.array(timem)
    obstime = Time(evt_time,format='mjd')#.to_value('isot')
    print(obstime)
    frame_hor = AltAz(obstime=obstime, location=locationAntares)
    local_sky=SkyCoord(az=azi,alt=alt,frame=frame_hor, unit=u.rad)
    print("after local sky")
    print("Now local to galactic")
    gal=local_sky.galactic#transform_to('galactic')
    print(gal.l.deg,gal.b.deg)
    return 0
    


if __name__ == "__main__":
    SPEED2coordinateconverter()
    
    
