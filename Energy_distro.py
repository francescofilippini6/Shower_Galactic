import numpy as np
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
from sklearn.metrics import confusion_matrix
import scipy

def reader(filename):
    print("Reading: ",filename)
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


def cut_dataframe_bdts(df):
    bdt_cut=0.12
    print("Slicing dataframe at: ", bdt_cut)
    selectedDF=df[df['BDT__cuts_1e2'] > bdt_cut]
    return selectedDF

def weight_astro_spectrum(df):
    spectral_index = 2.0
    spectral_norm = 6 * 10**-12 * 10**(3*spectral_index)  #Neutrino Energy in GeV
    #new_w3_weight=np.array(df['w2'])*spectral_norm*np.power(np.array(np.power(10,df['MCE'])),-spectral_index)
    my_spectral_index = 2.4
    my_spectral_norm = 2.1 * 10**-12 * 10**(3*my_spectral_index)  #Neutrino Energy in GeV
    #imporvements of cpu time of a factor 1000 minimum !!!!
    #new_w3_weight=np.array(df['w2'])*my_spectral_norm*np.power(np.array(np.power(10,df['MCE'])),-my_spectral_index)
    #error_norm=(4.9 *10**-12 * 10**(3*my_spectral_index))*np.power(np.array(np.power(10,df['MCE'])),-my_spectral_index)*np.array(df['w2']) #energy of nu  in GeV
    new_w3_weight=(np.array(df['w2'])*(4.8*10**-7))/(np.power(np.array(np.power(10,df['MCE'])),2.3))*(10**4)*(0.5)
    df['new_w3'] = np.array(new_w3_weight)
    #df['error_norm'] = np.array(error_norm)
    print("done")
    return df

def original_energy_distro(df,dfmuon):
    livetime=3012/365
    fig = plt.figure()
    ax=fig.add_subplot(111)
    #ax.hist(df['TantraEnergy'],bins=50,histtype='step',weights=np.array(df['new_w3'])*livetime,label='nu cosmic total')
    #print("first cut")
    df0=df[df['label']<0.1]
    print("first histo")
    ax.hist(df0['TantraEnergy'],bins=50,histtype='step',weights=np.array(df0['WeightAtmo']),label='nu atmo,not nu mu CC')
    print("second histo")
    ax.hist(df['TantraEnergy'],bins=50,histtype='step',weights=np.array(df['new_w3']),label= r'$nu cosmic, \gamma = -2.3$')
    print("second cut")
    df1=df[df['predicted_label']>0.8]
    print("third histo")
    ax.hist(df1['TantraEnergy'],bins=50,histtype='step',weights=np.array(df1['WeightAtmo']),label='nu atmo, numuCC')
    ax.hist(dfmuon['TantraEnergy'],bins=50,histtype='step',weights=np.array(dfmuon['WeightAtmo']),label='Atmo Muons')
    ax.set_title("Original Energy distribution, BDT>0.33")
    ax.set_xlabel('log10(E/GeV)')
    ax.set_ylabel(r'$\frac{dN}{dE}$')
    ax.set_yscale('log')
    ax.legend()
    plt.show()
    return fig

def original_energy_distroFEDE(df,dfmuon):
    livetime=3012/365
    fig = plt.figure()
    ax=fig.add_subplot(111)
    #ax.hist(df['TantraEnergy'],bins=50,histtype='step',weights=np.array(df['new_w3'])*livetime,label='nu cosmic total')
    #print("first cut")
    df0=df[df['interaction_type'].str.contains('numu')]
    print("first histo")
    ax.hist(df0['TantraEnergy'],bins=50,histtype='step',weights=np.array(df0['WeightAtmo']),label=r'$\nu^{\mu}_{atmo}$')
    print("second histo")
    ax.hist(df['TantraEnergy'],bins=50,histtype='step',weights=np.array(df['new_w3']),label= r'$nu cosmic, \gamma = -2.3$')
    print("second cut")
    df1=df[df['interaction_type'].str.contains('nue')]
    print("third histo")
    ax.hist(df1['TantraEnergy'],bins=50,histtype='step',weights=np.array(df1['WeightAtmo']),label=r'$\nu^e_{atmo}$')
    #ax.hist(dfmuon['TantraEnergy'],bins=50,histtype='step',weights=np.array(dfmuon['WeightAtmo'])*livetime,label='Atmo Muons')
    ax.set_title("Original Energy distribution, BDT>0.33")
    ax.set_xlabel('log10(E/GeV)')
    ax.set_ylabel(r'$\frac{dN}{dE}$')
    ax.set_yscale('log')
    #ax.set_ylim((10**-2,10**4))
    #ax.set_xlim((0.5,5.5))
    ax.legend()
    plt.show()
    return fig


def plot_reco_energy_distro_prediction0(df,dfmuon):
    livetime=3012/365
    fig = plt.figure()
    ax=fig.add_subplot(111)
    #ax.hist(df['TantraEnergy'],bins=50,histtype='step',weights=np.array(df['new_w3'])*livetime,label='nu cosmic total')
    print("first cut")
    df0=df[df['predicted_label']<0.1]
    print("first histo")
    ax.hist(df0['TantraEnergy'],bins=50,histtype='step',weights=np.array(df0['WeightAtmo']),label='nu atmo, prediction 0, evt number:'+ str(sum(df0['WeightAtmo'])))
    df1=df0[df0['label']<0.1]
    ax.hist(df1['TantraEnergy'],bins=50,histtype='step',weights=np.array(df1['WeightAtmo']),label='nu atmo, prediction 0, label 0, evt number:'+ str(sum(df1['WeightAtmo'])))
    df1=df0[df0['label']>0.8]
    ax.hist(df1['TantraEnergy'],bins=50,histtype='step',weights=np.array(df1['WeightAtmo']),label='nu atmo, prediction 0, label 1, evt number:'+ str(sum(df1['WeightAtmo'])))
    print("second histo")
    ax.hist(df0['TantraEnergy'],bins=50,histtype='step',weights=np.array(df0['new_w3']),label='nu cosmic, prediction 0, evt number:'+ str(sum(df0['new_w3'])))
    print("second cut")
    dfmuon=dfmuon[dfmuon['predicted_label']<0.1]
    ax.hist(dfmuon['TantraEnergy'],bins=50,histtype='step',weights=np.array(dfmuon['WeightAtmo']),label='Atmo Muons, prediction 0, evt number:'+ str(sum(dfmuon['WeightAtmo'])))
    ax.set_title('Predicted 0')
    ax.set_xlabel('log10(E/GeV)')
    ax.set_ylabel(r'$\frac{dN}{dE}$')
    ax.set_yscale('log')
    ax.legend()
    plt.show()
    return fig

def plot_reco_energy_distro_prediction1(df,dfmuon):
    livetime=3012/365
    fig = plt.figure()
    ax=fig.add_subplot(111)
    #ax.hist(df['TantraEnergy'],bins=50,histtype='step',weights=np.array(df['new_w3'])*livetime,label='nu cosmic total')
    print("first cut")
    df0=df[df['predicted_label']>0.8]
    ax.hist(df0['TantraEnergy'],bins=50,histtype='step',weights=np.array(df0['WeightAtmo'])*livetime,label='nu atmo, prediction 1, evt number:'+ str(sum(df0['WeightAtmo'])))
    df1=df0[df0['label']<0.1]
    ax.hist(df1['TantraEnergy'],bins=50,histtype='step',weights=np.array(df1['WeightAtmo'])*livetime,label='nu atmo, prediction 1, label 0, evt number:'+ str(sum(df1['WeightAtmo'])))
    df1=df0[df0['label']>0.8]
    ax.hist(df1['TantraEnergy'],bins=50,histtype='step',weights=np.array(df1['WeightAtmo'])*livetime,label='nu atmo, prediction 1, label 1, evt number:'+ str(sum(df1['WeightAtmo'])))
    print("second histo")
    ax.hist(df0['TantraEnergy'],bins=50,histtype='step',weights=np.array(df0['new_w3'])*livetime,label='nu cosmic, prediction 1, evt number:'+ str(sum(df0['new_w3'])))
    print("second cut")
    dfmuon=dfmuon[dfmuon['predicted_label']>0.8]
    ax.hist(dfmuon['TantraEnergy'],bins=50,histtype='step',weights=np.array(dfmuon['WeightAtmo']),label='Atmo Muons, prediction 0, evt number:'+ str(sum(dfmuon['WeightAtmo'])))
    ax.set_title('Predicted 1')
    ax.set_xlabel('log10(E/GeV)')
    ax.set_ylabel(r'$\frac{dN}{dE}$')
    ax.set_yscale('log')
    ax.legend()
    plt.show()
    return fig


def plot_reco_energy_distro_prediction0_checkdata(df,dfmuon,dfdata):
    livetime=3012/365
    fig = plt.figure()
    ax=fig.add_subplot(111)
    #ax.hist(df['TantraEnergy'],bins=50,histtype='step',weights=np.array(df['new_w3'])*livetime,label='nu cosmic total')
    print("first cut")
    df0=df[df['predicted_label']<0.1]
    print("first histo")
    ax.hist(df0['TantraEnergy'],bins=50,histtype='step',weights=np.array(df0['WeightAtmo']),label='nu atmo, prediction 0, evt number:'+ str(sum(df0['WeightAtmo'])))
    df1=df0[df0['label']<0.1]
    ax.hist(df1['TantraEnergy'],bins=50,histtype='step',weights=np.array(df1['WeightAtmo']),label='nu atmo, prediction 0, label 0, evt number:'+ str(sum(df1['WeightAtmo'])))
    df1=df0[df0['label']>0.8]
    ax.hist(df1['TantraEnergy'],bins=50,histtype='step',weights=np.array(df1['WeightAtmo']),label='nu atmo, prediction 0, label 1, evt number:'+ str(sum(df1['WeightAtmo'])))
    print("second histo")
    ax.hist(df0['TantraEnergy'],bins=50,histtype='step',weights=np.array(df0['new_w3']),label='nu cosmic, prediction 0, evt number:'+ str(sum(df0['new_w3'])))
    print("second cut")
    dfmuon=dfmuon[dfmuon['predicted_label']<0.1]
    ax.hist(dfmuon['TantraEnergy'],bins=50,histtype='step',weights=np.array(dfmuon['WeightAtmo']),label='Atmo Muons, prediction 0, evt number:'+ str(sum(dfmuon['WeightAtmo'])))


    dfdata=dfdata[dfdata['predicted_label']<0.5]
    data,bins=np.histogram(dfdata['TantraEnergy'],bins=50,weights=np.array(dfdata['WeightAtmo'])*10)
    ax.plot((bins[:-1]+bins[1:])/2,data,'ok')
    #ax.hist(dfdata['TantraEnergy'],bins=50,histtype='step',weights=np.array(dfdata['WeightAtmo'])*10,label='DATA'+ str(sum(dfdata['WeightAtmo'])))


    ax.set_title('Predicted 0')
    ax.set_xlabel('log10(E/GeV)')
    ax.set_ylabel(r'$\frac{dN}{dE}$')
    ax.set_yscale('log')
    ax.legend()
    plt.show()
    return fig



def confusion_matrix_weighted(df1):
    fig=plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('Total MC sample BDT>0.12')
    print("calculating confusion matrix")
    df_cm=confusion_matrix(df1['label'],df1['predicted_label'],sample_weight=df1['WeightAtmo'])
    df_cm=confusion_matrix(df1['label'],df1['predicted_label'],sample_weight=df1['new_w3'])
    sns.heatmap(df_cm, annot=True,cmap=plt.cm.Blues, fmt='g',ax=ax)
    print("saving figure")
    fig.savefig('Cosmic_weighted_2.3_BDT_0.12_confusion_matrix.png')
    plt.show()
    return fig

def angular_uncertanty(df):
    fig = plt.figure()
    df0=df[df['label']<0.1]
    ax=fig.add_subplot(221)
    #ax.plot(df0['MCZenith'],df0['TantraZenith'],'+r')
    ax.hist2d(df0['MCZenith'],df0['TantraZenith'], bins=50)
    ax.set_title('Shower-zenith')
    ax.set_xlabel('mc zenith')
    ax.set_ylabel('Tantra Zenith')
    #Print(df['TantraAngularEstimator'])
    #ax.hist(df['TantraAngularEstimator'],bins=50,histtype='step',label='Tantra Angular estimator distriution')
    #ax.set_title('Tantra Angular Estimator')
    #ax.set_xlabel('Angular Uncertanty')
    #ax.legend()
    df1=df[df['label']<0.8]
    ax2=fig.add_subplot(222)
    ax2.hist2d(df1['MCZenith'],df1['TantraZenith'], bins=50)
    ax2.set_title('Track-zenith')
    ax2.set_xlabel('mc zenith')
    ax2.set_ylabel('Tantra Zenith')
    
    ax3=fig.add_subplot(223)
    #ax.plot(df0['MCZenith'],df0['TantraZenith'],'+r')
    ax3.hist2d(df0['MCAzimuth'],df0['TantraAzimuth'], bins=50)
    ax3.set_title('Shower-Azimuth')
    ax3.set_xlabel('mc azimuth')
    ax3.set_ylabel('Tantra Azimuth')
    
    ax4=fig.add_subplot(224)
    #ax.plot(df0['MCZenith'],df0['TantraZenith'],'+r')
    ax4.hist2d(df1['MCAzimuth'],df1['TantraAzimuth'], bins=50)
    ax4.set_title('Track-Azimuth')
    ax4.set_xlabel('mc azimuth')
    ax4.set_ylabel('Tantra Azimuth')
    
    #ax.plot(df0['MCZenith'],df0['AAZenith'],'+b')
    #ax2.hist(df['TantraAngularEstimator'],bins=50,histtype='step',cumulative=True,label='Cumulative')
    #ax.legend()
    plt.show()
    return fig

if __name__ == "__main__":
    filename=sys.argv[1]
    muon=sys.argv[2]
    data=sys.argv[3]
    df=reader(filename)
    muon=reader(muon)
    data=reader(data)
    print("LENGTH:",len(df['TriggerT3']))
    df=df.drop(['TantraLines', 'TantraHits', 'Trigger3N', 'TriggerT3', 'Nrun',
       'TriggerCounter', 'FrameIndex', 'MCNGenEvent', 'MCNuType',
       'MCLeadingParticle', 'MCMultiplicity', 'Mestimator', 'TantraZenith',
       'TantraAzimuth', 'TantraAngularEstimator', 'TantraX', 'TantraY',
       'TantraZ', 'NOnTime', 'GridQuality', 'Lambda', 'AAZenith', 'AAAzimuth',
                 'Beta', 'TrackLength'],axis=1)
    #print(Counter(df['label']))
    #print(Counter(df['interaction_type']))
    muoncut=cut_dataframe_bdts(muon)
    df=cut_dataframe_bdts(df)
    datacut=cut_dataframe_bdts(data)
    print(df.keys())
    df1=weight_astro_spectrum(df)    
    
    print('AstroWeight sum:',sum(df1['new_w3']))
    print('AtmoWeight sum:',sum(df1['WeightAtmo']))
    #angular_uncertanty(df1)
    #original_energy_distro(df1,muoncut)
    #plot_reco_energy_distro_prediction1(df1,muon)
    #plot_reco_energy_distro_prediction0(df1,muon)
    plot_reco_energy_distro_prediction0_checkdata(df1,muoncut,datacut)
    
    
