import math as m
import numpy as np
import csv
import matplotlib.pyplot as plt
import astropy.units as u
import pandas as pd
import random as rand
import astropy.coordinates as coord
import astropy.units as u
from astropy.io import ascii
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, concatenate
from astropy import coordinates as coord
from astropy.coordinates.tests.utils import randomly_sample_sphere
from astropy.time import Time
from astropy import units as u
#from skyfield.api import load
#from skyfield.api import Topos
from collections import OrderedDict

import numpy as np
import pandas as pd

from km3astro.coord import local_event, local_frame, neutrino_to_source_direction
from matplotlib import pyplot as plt
import numpy as np
import healpy as hp

from km3astro.coord import (
    local_event,
    neutrino_to_source_direction,
    sun_local,
    get_location,
    convergence_angle,
    utm_zone,
    longitude_of_central_meridian,
)


x_bins = np.linspace(0,360,181)
y_bins = np.linspace(-90,90,91)
raa=[]
decl=[]
#signall=[]
#signalb=[]
for a  in range(100000):
    ra_deg=rand.uniform(0,360)
    decl_deg=rand.uniform(-61,8)
    #signall.append(rand.uniform(-40,40))
    #signalb.append(rand.uniform(-3,3))
    #if abs(ra_deg)<3 and abs(decl_deg)<40:
    #    continue
    raa.append(ra_deg)
    decl.append(decl_deg)
print(len(raa))
c = SkyCoord(raa, decl, frame="icrs", unit="deg")
#evt = local_event(azimuth, obstime, zenith, location="antares")
gal=c.galactic
gal_l_bkg_only=[]
gal_b_bkg_only=[]
signal_l=[]
lsignal_b=[]
#ciccio=[]
#ciccia=[]
#abs(gal.l[i])<40*u.deg
for i in range(len(gal.l)):
    ciccio=gal.l[i]< 40 *u.deg or gal.l[i]> 320 *u.deg
    if ciccio  and abs(gal.b[i])<3*u.deg:
        signal_l.append(gal.l[i])
        signal_b.append(gal.b[i])
        continue
    gal_l_bkg_only.append(gal.l[i])
    gal_b_bkg_only.append(gal.b[i])

bkg=SkyCoord(gal_l_bkg_only, gal_b_bkg_only, frame="galactic", unit="deg")
signal=SkyCoord(signal_l, signal_b, frame="galactic", unit="deg")
#signal=SkyCoord(signall, signalb, frame="galactic", unit="deg")
equatorial=bkg.transform_to("icrs")
galactic_ridge=signal.transform_to("icrs")
print(len(galactic_ridge.ra))
#--------------------------------------------------------------
# Setting of healpy resolution map and n. of pixels
#--------------------------------------------------------------
NSIDE = 64 # this sets the side resolution of each one of the 12 basic-macro pixels
print("Approximate resolution at NSIDE {} is {:.2} deg".format(NSIDE, hp.nside2resol(NSIDE, arcmin=True) / 60))
NPIX = hp.nside2npix(NSIDE)
print("number of pixels in the map",NPIX)

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
#--------------------------------------------------------------
#setting to zero all the map m
#--------------------------------------------------------------
m = np.zeros(hp.nside2npix(NSIDE))
hpx_map = cat2hpx(equatorial.ra.deg, equatorial.dec.deg, nside=64, radec=False)
hpx_map = cat2hpx(bkg.l.deg, bkg.b.deg, nside=32, radec=False)
#hpx_map = cat2hpx(eqsky.ra.deg,eqsky.dec.deg, nside=32, radec=False)
#hp.mollview(hpx_map, title="Mollview image RING")
hp.mollview(hpx_map, title="Antares bkg only - sky in galactic coord")
hp.graticule()
#plt.subplot(111, projection='aitoff')
#plt.grid(True)
#plt.plot(galactic_ridge.ra.wrap_at('180d').radian,galactic_ridge.dec.radian,'o',markersize=0.5, alpha=0.1, color='b')
plt.show()



