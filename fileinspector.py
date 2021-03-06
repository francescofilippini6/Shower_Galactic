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
import seaborn as sns
from matplotlib.gridspec import GridSpec
#from km3astro.coord import (
#    local_event,
#    neutrino_to_source_direction,
#    sun_local,
#    get_location,
#    convergence_angle,
#    utm_zone,
#    longitude_of_central_meridian,
#)

def reader(filename):
    df = pd.read_hdf(filename)
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
def expo(x,a,b):
    return a*np.power(x,b)
    
def plotBDT(listofdf,filename,option):
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
    summ=list(itertools.chain.from_iterable(sumhisto))
    summw=list(itertools.chain.from_iterable(sumweights))
    sumcosmicww=list(itertools.chain.from_iterable(sumcosmicweights))
    sumco=list(itertools.chain.from_iterable(sumcosmic))
    ax.hist(summ,bins=50,label='SUM',histtype=u'step',weights=summw)
    ax.hist(sumco,bins=50,label='cosmic',histtype=u'step',weights=sumcosmicww)
    plt.ylim((10**-2,10**6))
    plt.xlabel('BDT score')
    plt.legend()
    plt.show()
    return figu

def cut_dataframe_bdts(df):
    bdt_cut=0.33
    selectedDF=df[df['BDT__cuts_1e2'] > 0.33]
    return selectedDF
    

def plot_MC_RECO(df):
    Reco_zenith=df['AAZenith']
    Reco_azimuth=df['AAAzimuth']
    MCZen=df['MCZenith']
    MCAz=df['MCAzimuth']
    Reco_Energy=df['TantraEnergy']
    MC_Energy=df['MCE']
    print(Reco_zenith)
    print(MCZen)
    fig=plt.figure()
    ax1=fig.add_subplot(231)
    ax1.hist(Reco_zenith,bins=50,label='Reco_Zenith '+str(len(Reco_zenith)),alpha=0.5)
    ax1.hist(MCZen,bins=50,label='MC_Zenith '+str(len(MCZen)),alpha=0.5)
    ax1.set_title('RECO-MC Zenith')
    ax1.legend()
    ax2=fig.add_subplot(232)
    ax2.hist(Reco_azimuth,bins=50,label='Reco_Azimuth '+str(len(Reco_azimuth)),alpha=0.5)
    ax2.hist(MCAz,bins=50,label='MC_Azimuth '+str(len(MCAz)),alpha=0.5)
    ax2.set_title('RECO-MC Azimuth')
    ax2.legend()
    
    ax3=fig.add_subplot(233)
    ax3.hist(Reco_Energy,bins=50,label='Reco_Energy '+str(len(Reco_Energy)),alpha=0.5)
    ax3.hist(MC_Energy,bins=50,label='MC_Energy '+str(len(MC_Energy)),alpha=0.5)
    ax3.set_title('RECO-MC Energy')
    ax3.legend()

    ax3=fig.add_subplot(234)
    ax3.plot(Reco_Energy,MC_Energy)
    
    ax3=fig.add_subplot(235)
    ax3.plot(MCAz,Reco_azimuth)

    ax3=fig.add_subplot(236)
    ax3.plot(MCZen,Reco_zenith)
    
    plt.ylabel('entries')
    plt.show()
    return fig
    
def data_extractor(df,filename):
    name=filename.split('.')
    df.to_csv(name[0]+'.csv',index=False)
    return df

def weight_astro_spectrum(df):
    spectral_index = 2.0
    spectral_norm = 6 * 10**-12 * 10**(3*spectral_index)  #Neutrino Energy in GeV
    #new_w3_weight=np.array(df['w2'])*spectral_norm*np.power(np.array(np.power(10,df['MCE'])),-spectral_index)
    my_spectral_index = 2.5
    my_spectral_norm = 3.5 * 10**-12 * 10**(3*my_spectral_index)  #Neutrino Energy in GeV
    #imporvements of cpu time of a factor 1000 minimum !!!!
    new_w3_weight=np.array(df['w2'])*my_spectral_norm*np.power(np.array(np.power(10,df['MCE'])),-my_spectral_index)
    error_norm=(4.9 *10**-12 * 10**(3*my_spectral_index))*np.power(np.array(np.power(10,df['MCE'])),-my_spectral_index)*np.array(df['w2']) #energy of nu in GeV  
    df['new_w3'] = np.array(new_w3_weight)
    df['error_norm'] = np.array(error_norm)
    print("done")
    return df



def plot_reco_energy_distro(listofdf,filename):
    livetime=3012/365
    fig = plt.figure()
    gs=GridSpec(2,2)
    ax3=fig.add_subplot(gs[0,0]) # First row, first column
    ax2=fig.add_subplot(gs[0,1]) # First row, second column
    ax=fig.add_subplot(gs[1,:]) # Second row, span all columns

    cosmic_weight_old_w3=[]
    cosmic_weight_old_fede=[]
    cosmic_weight_new=[]
    cosmic_all=[]
    atmo_nu_e=[]
    atmo_nu_e_w=[]
    atmo_nu_mu=[]
    atmo_nu_mu_w=[]

    error_norm=[]
    
    for counter, df in enumerate(listofdf):
        if 'DATA' in filename[counter]:
            continue
        elif 'nue' in filename[counter]:
            atmo_nu_e.append(list(df['TantraEnergy']))
            atmo_nu_e_w.append(list(df['WeightAtmo']))
        
        elif 'numu' in filename[counter]:
            atmo_nu_mu.append(list(df['TantraEnergy']))
            atmo_nu_mu_w.append(list(df['WeightAtmo']))
            
        #atmoWeight = list(df['WeightAtmo'])
        #plotting the energy distribution of atmo events
        #ax.hist(df['TantraEnergy'],bins=50,histtype=u'step',weights=atmoWeight,label='E'+str(filename[counter].split('.')[0]))

        #collecting the energy values
        cosmic_all.append(list(df['TantraEnergy']))

        #collecting all the weights to produce the cosmic spectrum
        if 'new_w3' in df.keys():
            cosmic_weight_new.append(list(df['new_w3']))
            error_norm.append(list(df['error_norm']))
            

        cosmic_weight_old_fede.append(list(df['WeightAstro']))
        cosmic_weight_old_w3.append(np.array(df['w3']))
        
    
    #sum atmo plot
    sumnue_ene=list(itertools.chain.from_iterable(atmo_nu_e))
    sumnumu_ene=list(itertools.chain.from_iterable(atmo_nu_mu))
    sumnueatmo_w=list(itertools.chain.from_iterable(atmo_nu_e_w))
    
    
    
    sumnumuatmo_w=list(itertools.chain.from_iterable(atmo_nu_mu_w))
    #/len(sumnumu_ene)
    ax.set_yscale('log')
    #comsic flux plotting
    sumcomsic_ene=list(itertools.chain.from_iterable(cosmic_all))
    if 'new_w3' in df.keys():
        sumcosmicw3_new=list(itertools.chain.from_iterable(cosmic_weight_new))
        sumerrornorm_new=list(itertools.chain.from_iterable(error_norm))
        print("error norm",len(sumerrornorm_new))
        weigh=np.array(sumcosmicw3_new)
        print("weigh len",len(weigh))
        #print("Ciccio")
        histtn,binssn,_=ax.hist(sumcomsic_ene,bins=50,histtype=u'step',weights=np.array(weigh)*livetime,label='E-2.5 spectrum')

        ax.hist(sumcomsic_ene,bins=50,histtype=u'step',weights=np.array(sumerrornorm_new)*livetime,label='U-L for E-2.5',linestyle=('dashed'))

        #error derived by the norm flux uncertanty!
        dy=0#np.array(1.4 *10**-12)*np.power(np.power(10,np.array(binssn[1:])+np.array(binssn[:-1]))*10**-3,-2.5)
        #dy=np.array(sumerrornorm_new)*livetime
        print(dy)
        plt.errorbar((np.array(binssn[1:])+np.array(binssn[:-1]))/2, np.array(histtn), yerr=dy, fmt='.k', color='black',ecolor='lightgray', elinewidth=3, capsize=0)
        ax2.plot((np.array(binssn[1:])+np.array(binssn[:-1]))/2, np.array(histtn)*np.power(np.power(10,(np.array(binssn[1:])+np.array(binssn[:-1]))/2),2.5),'r+')
        ax2.set_xlabel('log10(E/GeV)')
        ax2.set_title('Distro my spectrum *E ^ 2.5')
        ax2.set_ylabel(r'$E^{2.5}*\frac{dN}{dE}$',rotation=90)
        ax2.set_ylim((-1000,10**4))
        
    sumcosmicw3_old_fede=list(itertools.chain.from_iterable(cosmic_weight_old_fede))
    #sumcosmicw3_old_w3=list(itertools.chain.from_iterable(cosmic_weight_old_w3))
    
    ax.hist(sumnumu_ene,bins=50,histtype=u'step',weights=np.array(sumnumuatmo_w),label='nu mu atmo')
    #/len(sumnumu_ene)
    ax.hist(sumnue_ene,bins=50,histtype=u'step',weights=np.array(sumnueatmo_w),label='nu e atmo')
    #/len(sumnue_ene)
    histt,binss,_=ax.hist(sumcomsic_ene,bins=50,histtype=u'step',weights=np.array(sumcosmicw3_old_fede),label='E-2.0 spectrum fede')
    #/len(sumcomsic_ene)
    ax.set_xlabel('log10(E/GeV)')
    ax.set_ylabel(r'$\frac{dN}{dE}$')
    ax.set_ylim((10**-2,10**4))
    ax.set_xlim((0.5,5.5))
    ax.legend()
    
    ax3.plot((np.array(binss[1:])+np.array(binss[:-1]))/2,np.array(histt)*np.power(np.power(10,(np.array(binss[1:])+np.array(binss[:-1]))/2),2.0),'b+')
    ax3.set_xlabel('log10(E/GeV)')
    ax3.set_title('Distro spectrum *E ^ 2.0')
    ax3.set_ylabel(r'$E^{2}*\frac{dN}{dE}$',rotation=90)
    ax3.set_ylim((-1000,10**4))
    plt.suptitle('Energy distributions')
    plt.show()
    
    print("FEDE COMSMIC EVENTS",np.sum(np.array(sumcosmicw3_old_fede)))
    print("MINE cosmic events",np.sum(np.array(weigh)*livetime))
    print("nu mu atmo",np.sum(np.array(sumnumuatmo_w)))
    print("nu e atmo",np.sum(np.array(sumnueatmo_w)))
    
    return fig
    

def histo_dataframe(df,filename):
    columns=df.keys()
    name=filename.split('.')
    fig=plt.figure(figsize=(20,10))
    fig.suptitle(name[0])
    for index,key in enumerate(columns):
        if any(sas == -99999 for sas in list(df[key])) or df[key].isnull:
            print('out')
            continue
        else:
            print('in')
            axx=fig.add_subplot(5,10,index+1)
            axx.hist(list(df[key]),bins=50,histtype=u'step')
    plt.show()
    return  fig

def histo_sn(df,filename):
    sns.pairplot(df,hue="species", diag_kind="hist")
    plt.show()

#def coordinateconverter(listofdataframe):
#    antares_lat=42.8  #42°48\' N
#    antares_lon=-6.17  #6°10' E ->  minus??
#    locationAntares = locationAntares= EarthLocation(lat=-antares_lat*u.deg , lon=(antares_lon+180)*u.deg, height= -12000*u.m)
#    for dataframe in listofdataframe:
#        evt_time=np.array(dataframe['DateMJD'])
#        print(evt_time)
#        #azi=np.degrees(np.array(dataframe['AAAzimuth']))
#        #conversion of the altitude to zeith!!!!
#        #alt=np.degrees(np.pi/2+np.array(dataframe['AAZenith']))
#        azi=np.degrees(np.array(dataframe['TantraZenith']))
#        alt=np.degrees(np.pi/2+np.array(dataframe['TantraAzimuth']))
#        gal_l=[]
#        gal_b=[]
#        print("alt (deg)",alt)
#        print("azi (deg)",azi)
#        for indexdf in range(len(alt)):
#            print("A",indexdf)
#            print(evt_time[indexdf])
#            obstime = Time(evt_time[indexdf],format='mjd').to_value('isot')
#            print("TIME",obstime)
#            frame_hor = AltAz(obstime=obstime, location=locationAntares)
#            print("After frame constructor")
#            print(azi[indexdf],alt[indexdf])
#            local_sky=SkyCoord(azi[indexdf],alt[indexdf],frame=frame_hor, unit=u.deg)
#            print("after local sky")
#            gal=local_sky.transform_to('galactic')
#            gal_l.append(gal.l)
#            gal_b.append(gal.b)
#            final=SkyCoord(l=gal_l,b=gal_b,frame='galactic', unit=u.deg)
#            dataframe['gal_coord']=final
#        return final

#onverter(listofdataframe):
#  #42°48\' N
#7  #6°10' E ->  minus??
#= locationAntares= EarthLocation(lat=-antares_lat*u.deg , lon=(antares_lon+180)*u.deg, height= -12000*u.m)
#)
# listofdataframe:
#
#
#ut the muon dataframe is empty, so keep attention and skip it !
#.empty:
#ciccio")
#
#array(dataframe['DateMJD'])
#me)
#ees(np.array(dataframe['TantraZenith']))
#ees(np.pi/2+np.array(dataframe['TantraAzimuth']))
#es(np.array(dataframe['AAAzimuth']))
#of the altitude to zeith!!!!
#es(np.pi/2+np.array(dataframe['AAZenith']))
#deg)",alt)
#deg)",azi)
#me(evt_time,format='mjd').to_value('isot')
#AltAz(obstime=obstime, location=locationAntares)
#ndexdf],alt[indexdf])
#yCoord(azi,alt,frame=frame_hor, unit=u.deg)
# local sky")
#y.transform_to('galactic')
# transformation")
#(gal.l)
#(gal.b)
#rd(l=gal_l,b=gal_b,frame='galactic', unit=u.deg)

def SPEED2coordinateconverter(listofdataframe):
    antares_lat=42.8  #42°48\' N
    antares_lon=-6.17  #6°10' E ->  minus??
    locationAntares = locationAntares= EarthLocation(lat=antares_lat*u.deg , lon= antares_lon *u.deg, height= -2500*u.m)
    print("entering")
    for dataframe in listofdataframe:
        #gal_l=[]
        #gal_b=[]
        #after the cut the muon dataframe is empty, so keep attention and skip it !
        if dataframe.empty:
            print ("ciccio")
            continue
        evt_time=np.array(dataframe['DateMJD'])
        print(evt_time)
        azi=np.array(dataframe['TantraAzimuth'])
        #it is called zenith but it is already the alt!!! Zenith (0,180), alt(-90,90)=90-zenith
        alt=np.array(dataframe['TantraZenith'])
        #azi=np.degrees(np.array(dataframe['AAAzimuth']))
        #alt=(np.array(dataframe['AAZenith']))
        print("alt (deg)",180*alt/np.pi)
        #print("azi (deg)",azi)
        obstime = Time(evt_time,format='mjd').to_value('isot')
        frame_hor = AltAz(obstime=obstime, location=locationAntares)
        #print(azi[indexdf],alt[indexdf])
        local_sky=SkyCoord(azi,alt,frame=frame_hor, unit=u.rad)
        print("after local sky")
        gal=local_sky.transform_to('galactic')
        print("after transformation")
        dataframe['gal_l'] = np.array(gal.l)
        dataframe['gal_b'] = np.array(gal.b)
        print(gal.l)
    return listofdataframe



def km3astro_coordinateconverter(listofdataframe):
    
    for dataframe in listofdataframe:
        gal_l=[]
        gal_b=[]
        #after the cut the muon dataframe is empty, so keep attention and skip it !
        if dataframe.empty:
            print ("ciccio")
            continue
        evt_time=np.array(dataframe['DateMJD'])
        print(len(evt_time))
        obstime = Time(evt_time,format='mjd').to_value('isot')
        #for t in obstime:
        #    print(t)
        #obstime = Time('2020-12-27T20:00')

        azimuth=np.array(dataframe['TantraAzimuth'])
        #con l'inserto di np.pi/2 non da più errore. Come mai?
        zenith=np.pi/2+np.array(dataframe['TantraZenith'])
        print(np.degrees(zenith))
        #for a in np.degrees(azimuth):
            #print(a)
        #for b in np.degrees(zenith):
            #print(b)
            #if b>90 or b<-90:
            #    print("MERDA")
        #azimuth, zenith = neutrino_to_source_direction(azimuth, zenith)
        evt = local_event(azimuth, obstime, zenith, location="antares")
        #print(evt)
        gal=evt.transform_to('galactic')
        #print(len(gal.l))
        dataframe['gal_l'] = np.array(gal.l)
        dataframe['gal_b'] = np.array(gal.b)

        #final=SkyCoord(l=gal.l,b=gal.b,frame='galactic', unit=u.deg)
    return listofdataframe

def galctic_ridge_event_selector(listofdataframe):
    for dataframe in listofdataframe:
        l_limit=40
        b_limit=3
        flag=[]
        for index in range(len(dataframe['gal_b'])):
            bkg=SkyCoord(dataframe['gal_l'][index]*180/np.pi, dataframe['gal_b'][index]*180/np.pi, frame="galactic", unit="deg")
            equatorial=bkg.transform_to("icrs")
            booleanL=dataframe['gal_l'][index]< 40 *u.deg or dataframe['gal_l'][index] > 320 *u.deg
            if booleanL and abs(dataframe['gal_b'])<b_limit*u.deg:
                flag.append('signal')
            elif equatorial.dec<8 and equatorial.dec>-61:
                flag.append('background')
            else:
                flag.append('outOR')
        dataframe['classifier']=np.array(flag)
    print(listofdataframe[0].keys())
    return listofdataframe

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
    option=sys.argv[2]
    cut=[]
    for a in filename:
        print('retrieving file: ',a)
        df=reader(a)
        histo_sn(df,a)
        if option=='cut':
            selection=cut_dataframe_bdts(df)
            cut.append(weight_astro_spectrum(selection))
            
        else:
            listofdf.append(weight_astro_spectrum(df))
        #if 'DATA' in a:
            #data_extractor(df,a)
            #histo_dataframe(df,a)
            #plotterskymap(coordinateconverter(selection))
         #   else:
         #       print('wait')
                #plotterskymap(coordinateconverter(df))
        #if 'anumuCC' in a:
            #plot_MC_RECO(df)
    print("CUT len",len(cut))
    print("--------GLOBAL----------")
    if option=='cut':
        #plot_reco_energy_distro(cut,filename)
        #GALACTIC=km3astro_coordinateconverter(cut)
        GALACTIC=SPEED2coordinateconverter(cut)
        galctic_ridge_event_selector(GALACTIC)
    else:
        plot_reco_energy_distro(listofdf,filename)
    #plotBDT(listofdf,filename,option)
    
 
