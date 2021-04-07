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
import seaborn as sns
import json
import h5py
import os
from sklearn.preprocessing import StandardScaler
from tensorflow.keras.models import model_from_json
from sklearn.metrics import confusion_matrix

#import tensorflow as tf
#print(tf.__version__)

def reader(filename):
    print("reading:", filename)
    df = pd.read_hdf(filename)
    return df

def cut_dataframe_bdts(df):
    #bdt_cut=0.33
    print("Bdt cut")
    bdt_cut=0.12
    selectedDF=df[df['BDT__cuts_1e2'] > bdt_cut]
    return selectedDF


def preprocessing(df1,df2):
    print("preprocessing dataset")
    aa = df1[['TantraLines', 'TantraHits', 'Mestimator', 'TantraZenith',
              'TantraAzimuth', 'TantraAngularEstimator', 'TantraX', 'TantraY',
              'TantraZ','Lambda','Beta', 'TrackLength','TantraEnergy','TantraRho','IntegralCharge','MeanCharge', 'StdCharge',
              'TriggerCounter','GridQuality','AAZenith', 'AAAzimuth','Trigger3N', 'TriggerT3','NOnTime']]
    features = df2[['TantraLines', 'TantraHits', 'Mestimator', 'TantraZenith',
              'TantraAzimuth', 'TantraAngularEstimator', 'TantraX', 'TantraY',
              'TantraZ','Lambda','Beta', 'TrackLength','TantraEnergy','TantraRho','IntegralCharge','MeanCharge', 'StdCharge',
              'TriggerCounter','GridQuality','AAZenith', 'AAAzimuth','Trigger3N', 'TriggerT3','NOnTime']]
    scaler = StandardScaler().fit(aa)
    #scaler = StandardScaler().fit(features)
    features = scaler.transform(features)
    #print("MEAN",features.mean(axis=0))
    #print("STD",features.std(axis=0))
    return features

def model_predicter(df):
    print("model loading...")
    json_file = open('model.json', 'r')
    loaded_model_json = json_file.read()
    loaded_model = model_from_json(loaded_model_json)
    # load weights into new model
    loaded_model.load_weights("best_model.hdf5")
    loaded_model.compile(optimizer='adam',loss='binary_crossentropy')
    print("model loaded")
    print("start prediction")
    y_pred=loaded_model.predict(df)
    return y_pred

def predicted_label_distribution(dfa,lista):
    notnumucc=[]
    numucc=[]
    labela=dfa.iloc[:,49:50].values
    labelb=np.concatenate(labela, axis=0)
    print("start the distribution")
    for a in range(len(labelb)):
        if labelb[a]<0.1:
            notnumucc.append(lista[a])
        elif labelb[a]>0.8:
            numucc.append(lista[a])
        else:
            print("Something wrong")
    print("finished splitting")
    histo1,bin_edges=np.histogram(notnumucc, bins=np.linspace(-0.5,1.05,30))#[0, 0.1, 0.2, 0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
    histo2=np.histogram(numucc, bins=bin_edges)
    print("histo1:", histo1)
    print("histo2:", histo2)
    return(histo1,histo2)
    #fig = plt.figure()
    #ax=fig.add_subplot(111)
    #ax.hist(notnumucc,bins=30,histtype='step',label='0 NN output distro')
    #ax.hist(numucc,bins=30,histtype='step',label='1 NN output distro')
    #plt.show()
    #return fig

if __name__ == "__main__":
    filename=sys.argv[1]
    scaler_distro=sys.argv[2]
    df=reader(filename)
    df1=reader(scaler_distro)
    if "MUON" in filename:
        df=cut_dataframe_bdts(df)
        print("MUON cut bdt")
    print("dataframe cut", len(df['TriggerT3']))
    predicted_labels=model_predicter(preprocessing(df1,df))
    #predicted_label_distribution(df,predicted_labels)
    y_pr=[]
    for a in predicted_labels:
        if a > 0.5:
            y_pr.append(1)
        else:
            y_pr.append(0)
     
    df['predicted_label']=predicted_labels
    #print("Classification",Counter(y_pr))
    print("storing result")
    df.to_hdf('Continuos_neutrino_prediction.h5', key='df', mode='w')

    #store = pd.HDFStore('HDF5_store_predicted.h5')
    #store.append('df',df)
    df_cm=confusion_matrix(np.ones(len(y_pr)),y_pr,sample_weight=df['WeightAtmo'])
    ax = plt.axes()
    ax.set_title(' sample')
    sns.heatmap(df_cm, annot=True,cmap=plt.cm.Blues, fmt='g',ax=ax)
    ax.set_ylabel('true')
    ax.set_xlabel('predicted')
    plt.show()
    
