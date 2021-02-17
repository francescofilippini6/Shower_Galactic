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
import healpy as hp
from scipy.optimize import curve_fit


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
def gaussian(x,B,a,b):
    return B*np.exp(-((x-a)**2/(2*b**2)))
    
def plotBDT(listofdf,namefile,option):
    figu, ax = plt.subplots()
    sumhisto=[]
    sumcosmic=[]
    sumweights=[]
    sumcosmicweights=[]
    for counter, df in enumerate(listofdf):
        if option=='cut':
            df=df[df['BDT__cuts_1e2'] > 0.33]
            if df.empty:
                continue
            print('cut')
        sumhisto.append(list(df['BDT__cuts_1e2']))
        sumweights.append(list(df['WeightAtmo']))

        #sumhisto=np.concatenate(sumhisto,np.array(df['BDT__cuts_1e2']),axis=None)
        aaa=np.array(df['WeightAtmo'])
        if 'DATA' in filename[counter]:
            weigh=[10]*len(df['BDT__cuts_1e2'])
            df['BDT__cuts_1e2'].plot.hist(ax=ax,bins=50,histtype=u'step',weights=weigh,label='DATA ('+str(len(df['BDT__cuts_1e2']))+' ponits)',color='k')
            #print("Number of Data points",len(df['BDT__cuts_1e2']))
        elif 'MUON' in filename[counter]:
            benfe,binning,_= ax.hist(df['BDT__cuts_1e2'],bins=50,histtype=u'step',weights=aaa,label=namefile[counter].split('.')[0])
            baricentri=np.array([0.5*(binning[i]+binning[i+1]) for i in range(len(binning)-1)])
            startRange=0.08
            baricentri = [i for i in baricentri if i >= startRange]
            popt,_=curve_fit(gaussian,baricentri,benfe[-len(baricentri):])
            x = np.linspace(startRange, max(df['BDT__cuts_1e2'])+0.2, 1000)
            ax.plot(x,gaussian(x,*popt))
        else:
            df['BDT__cuts_1e2'].plot.hist(ax=ax,bins=50, histtype=u'step',weights=aaa,label=namefile[counter].split('.')[0])
            sumcosmicweights.append(list(df['WeightAstro']))
            sumcosmic.append(list(df['BDT__cuts_1e2']))
        ax.set_yscale('log')
        #print(len(sumhisto[counter]))
    summ=list(itertools.chain.from_iterable(sumhisto))
    summw=list(itertools.chain.from_iterable(sumweights))
    sumcosmicww=list(itertools.chain.from_iterable(sumcosmicweights))
    sumco=list(itertools.chain.from_iterable(sumcosmic))
    #eee=np.array(listofdf[1]['WeightAtmo'])
    ax.hist(summ,bins=50,label='SUM',histtype=u'step',weights=summw)
    ax.hist(sumco,bins=50,label='cosmic',histtype=u'step',weights=sumcosmicww)
    #ax.hist(summ,bins=50,label='comsic',weights=)
    plt.ylim((10**-2,10**6))
    plt.xlabel('BDT score')
    plt.legend()
    plt.show()
    return figu

def plot_events_after_cuts(df):
    bdt_cut=0.33
    selectedDF=df[df['BDT__cuts_1e2'] > 0.33]
    print("selected",selectedDF)
    return selectedDF
    
def coordinateconverter(dataframe):
    antares_lat=42.8  #42°48\' N
    antares_lon=-6.17  #6°10' E ->  minus??
    locationAntares =locationAntares= EarthLocation(lat=-antares_lat*u.deg , lon=(antares_lon+180)*u.deg, height= -12000*u.m)
    evt_time=list(dataframe['DateMJD'])
    print(evt_time)
    azi=np.degrees(np.array(dataframe['AAAzimuth']))
    #conversion of the altitude to zeith!!!!
    alt=np.degrees(np.pi/2+np.array(dataframe['AAZenith']))
    gal_l=[]
    gal_b=[]
    print("alt (deg)",alt)
    print("azi (deg)",azi)
    for indexdf in range(len(alt)):
        print("A",indexdf)
        print(evt_time[indexdf])
        obstime = Time(evt_time[indexdf],format='mjd').to_value('isot')
        print("TIME",obstime)
        frame_hor = AltAz(obstime=obstime, location=locationAntares)
        print("After frame constructor")
        print(azi[indexdf],alt[indexdf])
        local_sky=SkyCoord(azi[indexdf],alt[indexdf],frame=frame_hor, unit=u.deg)
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
    hp.visufunc.mollview(cat2hpx(final.l.wrap_at('180d').deg, final.b.deg, nside=32, radec=False),hold=True,title='selected_data')
    hp.graticule()
    plt.show()
    return fig

#plotterskymap(coordinateconverter(df))

if __name__ == "__main__":
    filename=sys.argv[1].split(',')
    listofdf=[]
    for a in filename:
        print('retrieving file: ',a)
        df=reader(a)
        listofdf.append(df)
        #plotterskymap(coordinateconverter(plot_events_after_cuts(df)))
        if 'DATA' in a:
            print(len(df['DateMJD']))
            selection=plot_events_after_cuts(df)
            print(len(selection['DateMJD']))
            plotterskymap(coordinateconverter(selection))
            #plotterskymap(coordinateconverter(df))
    option=sys.argv[2]
    plotBDT(listofdf,filename,option)
    
 
