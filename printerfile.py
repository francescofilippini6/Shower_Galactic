import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math as m
import sys
import astropy.coordinates as coord
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, concatenate
from astropy.time import Time
import healpy as hp



def reader(filename):
    df = pd.read_hdf(filename)
    #print(df.keys())
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

def plotBDT(listofdf,namefile):
    figu, ax = plt.subplots()
    for counter, df in enumerate(listofdf):
        #print('atmo weight',np.sum(df['WeightAtmo']))
        #print('BDT',df['BDT__cuts_1e2'])
        #weightsum=[]
        aaa=np.array(df['WeightAtmo'])#*np.array(df['w2'])#*np.array(df['w3'])
        bbb=np.array(df['WeightAtmo'])*np.array(df['w2'])
        #for a in range(len(df['WeightAtmo'])):
        muon=[]
        if 'DATA' in filename[counter]:
            df['BDT__cuts_1e2'].plot.hist(ax=ax,bins=50,histtype=u'step',weights=aaa,label='DATA',color='k')
        elif 'MUON' in filename[counter]:
            df['BDT__cuts_1e2'].plot.hist(ax=ax,bins=50,histtype=u'step',weights=aaa,label=namefile[counter].split('.')[0])     
        else:
            df['BDT__cuts_1e2'].plot.hist(ax=ax,bins=50, histtype=u'step',weights=bbb,label=namefile[counter].split('.')[0])

        ax.set_yscale('log')
    plt.xlabel('BDT score')
    plt.legend()
    plt.show()
    return figu


def overalldistribution(df):
    fig, ax = plt.subplots()
    df.plot.hist(ax=ax,bins=550,histtype=u'step')
    plt.show()
    return fig

def coordinateconverter(df):
    antares_lat=42.8  #42°48\' N
    antares_lon=-6.17  #6°10' E ->  minus??
    locationAntares =locationAntares= EarthLocation(lat=-antares_lat*u.deg , lon=(antares_lon+180)*u.deg, height= -12000*u.m)
    evt_time=df['DateMJD']
    azi=np.degrees(np.array(df['AAAzimuth']))
    #conversion of the altitude to zeith!!!!
    alt=np.degrees(np.pi/2+np.array(df['AAZenith']))
    gal_l=[]
    gal_b=[]
    print("alt (deg)",alt)
    print("azi (deg)",azi)

    for a in range(10000):#range(len(alt)):
        print("A",a)
        obstime = Time(evt_time[a],format='mjd').to_value('isot')
        print("TIME",obstime)
        frame_hor = AltAz(obstime=obstime, location=locationAntares)
        print("After frame constructor")
        print(azi[a],alt[a])
        local_sky=SkyCoord(azi[a],alt[a],frame=frame_hor, unit=u.deg)
        print("after local sky")
        gal=local_sky.transform_to('galactic')
        gal_l.append(gal.l)
        gal_b.append(gal.b)

    final=SkyCoord(l=gal_l,b=gal_b,frame='galactic', unit=u.deg)
    return final


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

def plotterskymap(final):
    #---------------------------------------------
    #plotting part
    #---------------------------------------------
    #--------------------------------------------------------------
    # Setting of healpy resolution map and n. of pixels
    #--------------------------------------------------------------
    NSIDE = 32 # this sets the side resolution of each one of the 12 basic-macro pixels
    print("Approximate resolution at NSIDE {} is {:.2} deg".format(NSIDE, hp.nside2resol(NSIDE, arcmin=True) / 60))
    NPIX = hp.nside2npix(NSIDE)
    print("number of pixels in the map",NPIX)
    m = np.zeros(hp.nside2npix(NSIDE))
    fig=plt.figure()
    m2=fig.add_subplot(111)
    hp.visufunc.mollview(cat2hpx(final.l.wrap_at('180d').deg, final.b.deg, nside=32, radec=False),hold=True,title=str(sys.argv[1]))
    hp.graticule()
    plt.show()
    return fig

#plotterskymap(coordinateconverter(df))

if __name__ == "__main__":
    filename=sys.argv[1].split(',')
    for a in filename:
        df=reader(a)
        print(df['w2'])
        print(df['w3'])
        print(df['WeightAtmo'])
        overalldistribution(df)
