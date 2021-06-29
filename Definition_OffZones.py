import numpy as np
from matplotlib import gridspec
#from healpy.newvisufunc import projview, newprojplot
import random as rand
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
import healpy as hp

def converterMC():
    antares_northing = 4742381.9
    antares_easting = 268221.6
    antares_height = -2500  # m     (guessed, probably similar to orca)
    antares_utm_zone_number = 32
    antares_utm_zone_letter = "N"
    antares_utm_zone = "{num}{let}".format(num=antares_utm_zone_number, let=antares_utm_zone_letter)
    antares_latitude, antares_longitude = utm.to_latlon(antares_easting, antares_northing, antares_utm_zone_number, antares_utm_zone_letter)
    print(antares_latitude,antares_longitude)
    locationAntares = locationAntares= EarthLocation(lat=antares_latitude*u.deg , lon= antares_longitude*u.deg, height= antares_height*u.m)
    total_alt=[]
    total_az=[]
    offset_b=+2
    offset_l=0
    for i in range(9):   #for i in range(10):
        shift=4.5/24+i*(1-4.5/24-7.2/24)/9
        obstime = Time(shift,format='mjd')
        b=[]
        l=[]
        for a in range(10000):
            aaa=rand.uniform(-3+offset_b,3+offset_b)
            ccc=rand.uniform(-40+offset_l,40+offset_l)
            if ccc<0:
                ccc+=360
            b.append(aaa)
            l.append(ccc)
        aa=SkyCoord(b=b,l=l,frame='galactic', unit=u.deg)
        #obb= Time(58965,format='mjd')
        frame_hor = AltAz(obstime=obstime,location=locationAntares)
        cc=aa.transform_to(frame_hor)
        #frame_horrr = AltAz(obstime=obstime,location=locationAntares)
        #kkk=cc.transform_to(frame_horrr)
        total_alt.append(cc.alt)
        total_az.append(cc.az)
    oo_time=Time(55500,format='mjd')
    ootime = Time(oo_time,format='mjd')
    frame_single = AltAz(obstime=ootime, location=locationAntares)
    local_sky=SkyCoord(az=total_az,alt=total_alt,frame=frame_single, unit=u.deg)
    ss=local_sky.transform_to('galactic')
    
    
    #for bb in range(10):
    #    print(ss[bb].dec.deg, ss[bb].ra.deg)
    #    print(max(ss[bb].dec.deg),max(ss[bb].ra.deg))
    #    print(min(ss[bb].dec.deg),min(ss[bb].ra.deg))
    return ss
    

def cat2hpx(lon, lat, nside, radec=True):
    """
    Convert a catalogue to a HEALPix map of number counts per resolution
    element.

    Parameters
    ----------
    lon, lat : (ndarray, ndarray)
        Coordinates of the sources in degree. If radec=True, assume input is in the icrs
        coordinate system. Otherwise assume input is glon, glat

    nside : int
        HEALPix nside of the target map

    radec : bool
        Switch between R.A./Dec and glon/glat as input coordinate system.
    Return
    ------
    hpx_map : ndarray
        HEALPix map of the catalogue number counts in Galactic coordinates
    
    """
    npix = hp.nside2npix(nside)
    
    if radec:
        eq = SkyCoord(lon, lat, 'icrs', unit='deg')
        l, b = eq.galactic.l.value, eq.galactic.b.value
    else:
        l, b = lon, lat
    # conver to theta, phi
    #theta conversion to port the b [-90,90] to the colatitude range [0,180]
    theta = np.radians(90. - b)
    phi = np.radians(l)
    # convert to HEALPix indices
    indices = hp.ang2pix(nside, theta, phi)
    
    idx, counts = np.unique(indices, return_counts=True)
    
    # fill the fullsky map
    hpx_map = np.zeros(npix, dtype=int)
    hpx_map[idx] = counts
    
    return hpx_map




if __name__ == "__main__":
    #galactic_ridge=converterMC()
    NSIDE = 64 # this sets the side resolution of each one of the 12 basic-macro pixels
    print("Approximate resolution at NSIDE {} is {:.2} deg".format(NSIDE, hp.nside2resol(NSIDE, arcmin=True) / 60))
    NPIX = hp.nside2npix(NSIDE)
    print("number of pixels in the map",NPIX)    
    #--------------------------------------------------------------
    #setting to zero all the map m
    #--------------------------------------------------------------
    m = np.zeros(hp.nside2npix(NSIDE))
    f=plt.figure()
    h1=f.add_subplot(211)
    hpx_map = cat2hpx(-360+galactic_ridge.l.deg,galactic_ridge.b.deg, nside=32, radec=False)
    hp.mollview(hpx_map, title="Antares all tracks, $\Lambda>5$,$\beta<1$",hold=True)
    hp.graticule()

    hh=f.add_subplot(212)
    hh.plot(galactic_ridge.l.deg,galactic_ridge.b.deg,'.',markersize=2,color='g')
    hh.axvline(x=220,color='k',linestyle='--')
    hh.axvline(x=330,color='k',linestyle='--')
    plt.show()
    
    

