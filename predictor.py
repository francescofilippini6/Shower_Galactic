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
    print("MEAN",features.mean(axis=0))
    print("STD",features.std(axis=0))
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


if __name__ == "__main__":
    filename=sys.argv[1]
    scaler_distro=sys.argv[2]
    df=reader(filename)
    df1=reader(scaler_distro)
    #df1=cut_dataframe_bdts(df)
    #df=df[df['interaction_type']=='numuCC']
    print("dataframe cut", len(df['TriggerT3']))
    predicted_labels=model_predicter(preprocessing(df1,df))
    y_pr=[]
    counter0=0
    counter1=0
    for a in predicted_labels:
        if a > 0.5:
            y_pr.append(1)
            counter1+=1
        else:
            y_pr.append(0)
            counter0+=1
    df['predicted_label']=y_pr
    print(df1.keys())
    df.to_hdf('Prediction_MC_all_data_BDT_CUT_0.12.h5', key='df', mode='w')
    df_cm=confusion_matrix(df['label'],y_pr)
    ax = plt.axes()
    ax.set_title('Total MC sample BDT>0.12')
    sns.heatmap(df_cm, annot=True,cmap=plt.cm.Blues, fmt='g',ax=ax)
    ax.set_ylabel('true')
    ax.set_xlabel('predicted')
    plt.show()
    #print("accuracy: ",counter0/(counter1+counter0))
