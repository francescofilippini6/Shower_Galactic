import numpy as np
from matplotlib import gridspec
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
from astropy.coordinates import FK5
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
from collections import Counter

zoneAmplitudeb=3
zoneAmplitudel=40
noff=9

def zonesSelection(dfa,lat,lon):
    print("cut on b latitude")
    dfb=dfa[np.absolute(dfa[lat])<zoneAmplitudeb]
    print("cut on l longitude")
    dfc=dfb[(dfb[lon]<zoneAmplitudel) | (dfb[lon]>360-zoneAmplitudel)]
    return dfc

def converter(df):
    antares_northing = 4742381.9
    antares_easting = 268221.6
    antares_height = -2500  # m     (guessed, probably similar to orca)
    antares_utm_zone_number = 32
    antares_utm_zone_letter = "N"
    antares_utm_zone = "{num}{let}".format(num=antares_utm_zone_number, let=antares_utm_zone_letter)
    antares_latitude, antares_longitude = utm.to_latlon(antares_easting, antares_northing, antares_utm_zone_number, antares_utm_zone_letter)
    print(antares_latitude,antares_longitude)
    locationAntares = locationAntares= EarthLocation(lat=antares_latitude*u.deg , lon= antares_longitude*u.deg, height= antares_height*u.m)
    azor=np.array(df['TantraAzimuth'])
    
    altor=np.array(np.pi/2-np.arccos(df['TantraZenith']))

    Timemjd=np.array(df['DateMJD'])
    obstime = Time(Timemjd,format='mjd')
    #-------- frame to be used icrs or FK5 (= j2000)
    frame_original = AltAz(obstime=obstime,location=locationAntares)
    c = SkyCoord(alt=altor, az=azor,frame=frame_original,unit=u.rad)
    print("Now local to galactic")
    gal=c.galactic
    df['gal_l']=gal.l.deg
    df['gal_b']=gal.b.deg
    return df

def checker_off_zone(df):
    antares_northing = 4742381.9
    antares_easting = 268221.6
    antares_height = -2500  # m     (guessed, probably similar to orca)
    antares_utm_zone_number = 32
    antares_utm_zone_letter = "N"
    antares_utm_zone = "{num}{let}".format(num=antares_utm_zone_number, let=antares_utm_zone_letter)
    antares_latitude, antares_longitude = utm.to_latlon(antares_easting, antares_northing, antares_utm_zone_number, antares_utm_zone_letter)
    print(antares_latitude,antares_longitude)
    locationAntares = locationAntares= EarthLocation(lat=antares_latitude*u.deg , lon=antares_longitude*u.deg, height= antares_height*u.m)
    #--------------------------------------------------
    azoff=np.array(df['TantraAzimuth'])
    altoff=np.array(np.pi/2-np.arccos(df['TantraZenith']))
    Timemjd=np.array(df['DateMJD'])
    #timecommon=Time(Timemjd,format='mjd')
    #frame_common = AltAz(obstime=timecommon,location=locationAntares)
    #ca = SkyCoord(az=azoff, alt=altoff,frame=frame_common,unit=u.rad)
    
    for i in range(noff+1):  
        shifttime=4.5/24. + i*(1 - 4.5/24 - 7.2/24 )/noff 
        shift=Timemjd+shifttime
        oobstime = Time(shift,format='mjd')
        frame_hor = AltAz(obstime=oobstime,location=locationAntares)
        cccc=SkyCoord(alt=altoff, az=azoff,frame=frame_hor,unit=u.rad)
        #print("second",cccc)
        gallo=cccc.galactic
        off=[]
        off=np.zeros(len(gallo.l))
        for coor in range(len(gallo.l)):
            if np.absolute(float(gallo[coor].b.deg))<zoneAmplitudeb:
                if float(gallo[coor].l.deg)<zoneAmplitudel or float(gallo[coor].l.deg)>(360-zoneAmplitudel):
                    off[coor]=1
                else:
                    continue
            else:
                continue
        df['off'+str(i)]=off
        df['gal_b'+str(i)]=gallo.b.deg
        df['gal_l'+str(i)]=gallo.l.deg
    return df

def energy_distro(df):
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.hist(10**df['log10rho'],histtype='step',bins=50,label='tracks PS-2020')
    ax.set_xlabel(r'log$_{10}\rho$')
    ax.legend()
    #ax.title('ON Zone distribution')
    print('MAX',max(df['log10rho']))
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.show()
    return fig


def plot_ONZONE(df,bins):
    print("ON Zone")
    #---------------------------------------------
    # cut on Tantra directions to be in ON
    #---------------------------------------------

    df=zonesSelection(df,'gal_b','gal_l')
    print("ON histo")
    hist0, _ = np.histogram(df['TantraEnergy'], bins=bins)
   
    #--------------------------------------------------
    #  also MC zenith and azimuth must be in the ON region
    #--------------------------------------------------
    
    print("ON events:", sum(hist0))
    return (hist0)

def plot_OFFZONE(df):
    print("OFF_Zones_selection")
    df0=zonesSelection(df,'gal_b0','gal_l0')
    df1=zonesSelection(df,'gal_b1','gal_l1')
    df2=zonesSelection(df,'gal_b2','gal_l2')
    df3=zonesSelection(df,'gal_b3','gal_l3')
    df4=zonesSelection(df,'gal_b4','gal_l4')
    df5=zonesSelection(df,'gal_b5','gal_l5')
    df6=zonesSelection(df,'gal_b6','gal_l6')
    df7=zonesSelection(df,'gal_b7','gal_l7')
    df8=zonesSelection(df,'gal_b8','gal_l8')
    #df9=zonesSelection(df,'gal_b9','gal_l9')
    #df10=zonesSelection(df,'gal_b10','gal_l10')
    #df11=zonesSelection(df,'gal_b11','gal_l11')
    print("OFF histo")

    
    hist0, bins = np.histogram(df0['TantraEnergy'], bins=50)
    hist1, _ = np.histogram(df1['TantraEnergy'], bins=bins)
    hist2, _ = np.histogram(df2['TantraEnergy'], bins=bins)
    hist3, _ = np.histogram(df3['TantraEnergy'], bins=bins)
    hist4, _ = np.histogram(df4['TantraEnergy'], bins=bins)
    hist5, _ = np.histogram(df5['TantraEnergy'], bins=bins)
    hist6, _ = np.histogram(df6['TantraEnergy'], bins=bins)
    hist7, _ = np.histogram(df7['TantraEnergy'], bins=bins)
    hist8, _ = np.histogram(df8['TantraEnergy'], bins=bins)
    #hist9, _ = np.histogram(df9[''TantraEnergy''], bins=bins)
    #hist10, _ = np.histogram(df10[''TantraEnergy''], bins=bins)
    #hist11, _ = np.histogram(df11[''TantraEnergy''], bins=bins)

    histsum=(hist0+hist1+hist2+hist3+hist4+hist5+hist6+hist7+hist8)/noff
    print("mean OFF events:",sum(histsum))
    #print(bins)
    return (histsum,bins)


def ONandOFF(df):
    offhisto,binning=plot_OFFZONE(df)
    onhisto=plot_ONZONE(df,binning)
    center = (binning[:-1] + binning[1:]) / 2
    cumulativeon=[]
    cumulativeoff=[]
    
    for i in range(len(binning)-1):
        cumulativeon.append(sum(onhisto[i:]))
        cumulativeoff.append(sum(offhisto[i:]))
    
        
    #saver.to_csv('Energy_Mrf.csv')
    #------------------------------------------------------
    # relative discrepancy
    #------------------------------------------------------
    #discrepancy=(np.absolute(onhisto-offhisto)/offhisto)*100
    
    #------------------------------------------------------
    # ratio discrepancy
    #------------------------------------------------------
    #discrepancy=onhisto/offhisto
    #------------------------------------------------------
    # sigma discrepancy
    #------------------------------------------------------
    
    discrepancy=[]
    for binc in range(len(onhisto)):
        tailpos=poisson.cdf(onhisto[binc],offhisto[binc]) #cdf=cumulative distribution function
        pvalue=1-tailpos
        discrepancy.append(st.norm.ppf(1-pvalue))

    discrepancy1=[]
    for bind in range(len(onhisto)):
        LieMa=np.sqrt(2)*np.sqrt(onhisto[bind]*np.log(2*onhisto[bind]/(onhisto[bind]+offhisto[bind]))+offhisto[bind]*np.log(2*offhisto[bind]/(onhisto[bind]+offhisto[bind])))
        discrepancy1.append(LieMa)
    #
    #discrepancy2=[]
    #for bine in range(len(onhisto)):
    #    onoffarticle=sc.betainc(0.5,onhisto[bine],1+offhisto[bine])/sc.beta(onhisto[bine],1+offhisto[bine])
    #    discrepancy2.append(onoffarticle)
    #print(discrepancy2)
    #print(discrepancy,len(discrepancy))

    fig = plt.figure()
    fig.suptitle('Energy distro On/OFF, l:'+str(zoneAmplitudel)+',b:'+str(zoneAmplitudeb), fontsize=16)
    gs = fig.add_gridspec(nrows=4, ncols=2,  hspace=0)#gs=GridSpec(4,2)
    ax=fig.add_subplot(gs[:-1,0])
    ax.plot(center,offhisto,'+r',label='offzone')
    ax.plot(center,onhisto,'+b',label='onzone')
    ax.set_xlabel(r'log$_{10}(TantraE/GeV)$')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.legend()
    ax1=fig.add_subplot(gs[:,1])
    ax1.plot(center,cumulativeoff,'--r',label='offzone')
    ax1.plot(center,cumulativeon,'--b',label='onzone')
    ax1.set_xlabel(r'log$_{10}(TantraEnergy/GeV)$')
    ax1.set_title('Cumulative')
    ax1.set_ylabel(r'$\frac{dN}{dE}$')
    ax1.set_yscale('log')
    ax1.legend()
    ax2=fig.add_subplot(gs[-1:,0])
    ax2.axhline(y=2, color='g', linestyle='-',label=r'2 $\sigma$')
    ax2.axhline(y=3, color='r', linestyle='-',label=r'3 $\sigma$')
    ax2.axhline(y=5, color='b', linestyle='-',label=r'5 $\sigma$')
    #ax2.plot(center,discrepancy,'+k',label=r'poisson ($\mu_s$=0)')
    ax2.plot(center,discrepancy1,'xb',label='Li and Ma')
    #ax2.plot(center,discrepancy2,'.r',label='Cousins Linnemann Tucker')
    ax2.legend(fontsize=7)
    ax2.set_ylabel(r'$\sigma$')
    ax2.set_xlabel(r'log$_{10}(TantraEnergy/GeV)$')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    plt.show()
    return (offhisto,onhisto)




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


def plotter(array_l,array_b,coordinate_frame):
    galactic_ridge=SkyCoord(l=array_l,b=array_b,frame=coordinate_frame, unit=u.deg)
    NSIDE = 64 # this sets the side resolution of each one of the 12 basic-macro pixels
    print("Approximate resolution at NSIDE {} is {:.2} deg".format(NSIDE, hp.nside2resol(NSIDE, arcmin=True) / 60))
    NPIX = hp.nside2npix(NSIDE)
    print("number of pixels in the map",NPIX)
    #setting to zero all the map m
    #--------------------------------------------------------------
    m = np.zeros(hp.nside2npix(NSIDE))
    f=plt.figure()
    h1=f.add_subplot(211)
    hpx_map = cat2hpx(galactic_ridge.l.deg,galactic_ridge.b.deg, nside=32, radec=False)
    hp.mollview(hpx_map, title=r"Antares OFF tracks, $\Lambda>5$,$\beta<1$",hold=True)
    hp.graticule()

    hh=f.add_subplot(212)
    hh.plot(galactic_ridge.l.deg,galactic_ridge.b.deg,'.',markersize=2,color='g')
    #hh.axvline(x=-360+220,color='k',linestyle='--')
    #hh.axvline(x=-360+330,color='k',linestyle='--')
    plt.show()
    return 0

if __name__ == "__main__":
    dfaa=pd.read_hdf('DataShowerPrediction.h5')
    dfaa=dfaa[dfaa['BDT__cuts_1e2']>0.1]
    dfaa=dfaa[dfaa['predicted_dropout']<0.3]
    #dfaa=dfaa[dfaa['TantraEnergy']>3.582]
    
    print('ZENITH',dfaa['TantraZenith'])
    print(min(dfaa['TantraZenith']),max(dfaa['TantraZenith']))
    
    dfconverted=converter(dfaa)
    dfbb=checker_off_zone(dfaa)
    dfbb.to_csv('Shower_2020.csv',mode='w')
    
    dfdd = pd.read_csv('Shower_2020.csv')
    print(dfdd.keys())
    lista=[]
    
    for a in range(len(dfdd['TantraZenith'])):
        aaaaa=dfdd.iloc[a]
        aaa=list(aaaaa[['off0','off1','off2','off3','off4','off5','off6','off7','off8']])
        counter=Counter(aaa)
        if counter[1]>1:
            lista.append(a)
        print('Len lista',len(lista))
    print('LISTA',lista)
    #modDfObj = dfdd.drop(lista)
    #print(modDfObj)
    ONandOFF(dfdd)    
    off=['off0','off1','off2','off3','off4','off5','off6','off7','off8']
    lv=[]
    bv=[]
    for b in off:
        for a in range(len(dfdd['TantraZenith'])):
            if dfdd[b].iloc[a] >0.9:
                lv.append(dfdd['gal_l'].iloc[a])
                bv.append(dfdd['gal_b'].iloc[a])
    print(len(lv))
    plotter(lv,bv,'galactic')
