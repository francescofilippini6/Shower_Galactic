# -*- coding: utf-8 -*-
"""MLP_binary_classificator.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1XUqQmNb0FIXCjll4YgkVBwG-cv_3EfQC
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json
import h5py
import os
from keras.models import Sequential,load_model
from keras.layers import Activation, Dropout, Flatten, Dense
from keras.models import Model
from keras.callbacks import ModelCheckpoint,Callback,EarlyStopping
from keras.optimizers import Adam,Adadelta, RMSprop
from keras.utils import plot_model, to_categorical
from keras import backend as K
from keras.utils.generic_utils import get_custom_objects
#from sklearn.preprocessing import Imputer
from keras.callbacks import LearningRateScheduler
from sklearn.preprocessing import StandardScaler #distributing the input values according to a norm distribution (men=0 std err=1) Scaling data
from sklearn.preprocessing import MinMaxScaler #normalize data x-min/max-min
from keras.callbacks import ReduceLROnPlateau # after n epochs in which the metric monitored do not improve the learning rate is reduced, to fine tune
#kernel_initializer="he_normal, init='uniform'.  he_uniform?
from sklearn.model_selection import train_test_split

def exp_decay(epoch):
    decay_rate = learning_rate / (1+epoch)
    lrate = learning_rate * np.exp(-decay_rate*epoch)
    return lrate
  
lr_rate = LearningRateScheduler(exp_decay)

checkpoint = ModelCheckpoint("best_model.hdf5", monitor='val_accuracy', verbose=1,save_best_only=True, mode='auto')
reduce_lr = ReduceLROnPlateau(monitor='val_accuracy', factor=0.2,patience=8, min_lr=0.001)
callbacks_list = [checkpoint,reduce_lr]

df = pd.read_hdf('/content/drive/MyDrive/dataset_BDT_CUT_0.12.h5')
#labels=df['label']
#mce=df['MCE']
labels=df.iloc[:,24:26]
features=df.iloc[:,0:24]
#features=features.drop(['AAZenith', 'AAAzimuth','Trigger3N', 'TriggerT3','IntegralCharge','MeanCharge', 'StdCharge','TriggerCounter'],axis=1) 
#adding MCE

print(features.keys())
print(labels.keys())
#replacing the null values with the mean of that column
#imputer = Imputer()
#X = imputer.fit_transform(X)
# Define the scaler 
scaler = StandardScaler().fit(features)
#scaler = MinMaxScaler().fit(features)
features = scaler.transform(features)

X_train, X_test, y_train, y_test = train_test_split(features, labels, test_size=0.2, random_state=42)


y_train=y_train['label']
print(y_test.keys())

#-------------------------------------------------
# cell for standard training
#-------------------------------------------------


model = Sequential()
model.add(Dense(64,activation="relu",kernel_initializer="he_normal", input_shape=(24,)))
#model.add(Dropout(0.35))
model.add(Dense(32,kernel_initializer="he_normal",activation="relu"))
#model.add(Dropout(0.35))
#model.add(Dense(,kernel_initializer="he_normal",activation="relu"))
model.add(Dense(32,kernel_initializer="he_normal",activation="relu"))
model.add(Dense(16,kernel_initializer="he_normal",activation="relu"))
model.add(Dense(1,activation='sigmoid'))
    


#adam2 = Adam(lr=learning_rate, momentum=momentum, decay=decay_rate, nesterov=False)
#adam2=Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay=0.0)
#adagrad, adadelta, RMSprop
model.compile(optimizer='adam',
              loss='binary_crossentropy',
              metrics=['accuracy'])
batch_size=32
nb_epoch=1
history = model.fit(X_train, y_train,             
          batch_size=batch_size,
          epochs=nb_epoch,
          validation_split=0.2,
          shuffle=True,
          verbose=1,
          callbacks=callbacks_list)


print(model.summary())

plt.plot(history.history['loss'])
print(history.history['loss'])
print(history.history['val_loss'])
plt.plot(history.history['val_loss'])
plt.title('model loss')
plt.ylabel('loss')
plt.xlabel('epoch')
plt.legend(['train', 'validation'], loc='upper left')
plt.show()

plt.plot(history.history['accuracy'])
print(history.history['accuracy'])
print(history.history['val_accuracy'])
plt.plot(history.history['val_accuracy'])
plt.title('model accuracy')
plt.ylabel('accuracy')
plt.xlabel('epoch')
plt.legend(['train', 'validation'], loc='upper left')
plt.show()

!pip install eli5

#-------------------------------------------------
# cell for feture importance extraction
#-------------------------------------------------

import eli5
from eli5.sklearn import PermutationImportance
from keras.wrappers.scikit_learn import KerasClassifier, KerasRegressor
from sklearn.utils.class_weight import compute_class_weight

def base_model():
    model = Sequential()
    model.add(Dense(64,activation="relu",kernel_initializer="he_normal", input_shape=(24,)))
    #model.add(Dropout(0.35))
    model.add(Dense(32,kernel_initializer="he_normal",activation="relu"))
    #model.add(Dropout(0.35))
    #model.add(Dense(,kernel_initializer="he_normal",activation="relu"))
    model.add(Dense(32,kernel_initializer="he_normal",activation="relu"))
    model.add(Dense(16,kernel_initializer="he_normal",activation="relu"))
    model.add(Dense(1,activation='sigmoid'))
    
    model.compile(optimizer='adam',
              loss='binary_crossentropy',
              metrics=['accuracy'])
    return model


#sk_params=[]

batch_size=32
nb_epoch=20
#my_model = KerasClassifier(build_fn=base_model)    
classifier = KerasClassifier(build_fn = base_model,validation_split=0.2,batch_size=batch_size,shuffle=True,epochs=nb_epoch,verbose=1,callbacks=callbacks_list)
classifier.fit(X_train, y_train)

perm = PermutationImportance(classifier, random_state=1,n_iter=1).fit(X_test, y_test)
eli5.show_weights(perm, feature_names = ['TantraLines', 'TantraHits', 'Mestimator', 'TantraZenith',
       'TantraAzimuth', 'TantraAngularEstimator', 'TantraX', 'TantraY',
       'TantraZ', 'Lambda', 'Beta', 'TrackLength', 'TantraEnergy', 'TantraRho',
       'IntegralCharge', 'MeanCharge', 'StdCharge', 'TriggerCounter',
       'GridQuality', 'AAZenith', 'AAAzimuth', 'Trigger3N', 'TriggerT3',
       'NOnTime'])

perm_train_feat_imp_df = pd.DataFrame({'val': perm.results_[0],
                                      'lab':['TantraLines', 'TantraHits', 'Mestimator', 'TantraZenith',
       'TantraAzimuth', 'TantraAngularEstimator', 'TantraX', 'TantraY',
       'TantraZ', 'Lambda', 'Beta', 'TrackLength', 'TantraEnergy', 'TantraRho',
       'IntegralCharge', 'MeanCharge', 'StdCharge', 'TriggerCounter',
       'GridQuality', 'AAZenith', 'AAAzimuth', 'Trigger3N', 'TriggerT3','NOnTime'] } )
perm_train_feat_imp_df.plot.barh(x='lab', y='val')

!pip install shap

import shap
background = X_train[np.random.choice(X_train.shape[0], 400, replace=False)]
explainer = shap.KernelExplainer(model, background)
shap_values = explainer.shap_values(X_test[0:20])

shap.initjs()
a=shap.force_plot(explainer.expected_value[0], 
                  shap_values[0][0,:], 
                  X_test[0,:])
display(a)

#shap_values50 = explainer.shap_values(X_test[280:330], nsamples=500)
shap.initjs()
shap.force_plot(explainer.expected_value, shap_values[0], X_test[10:])

shap.summary_plot(shap_values, X_test[0:10,:])

from keras.models import model_from_json
model_json = model.to_json()
with open("model.json", "w") as json_file:
    json_file.write(model_json)
json_file = open('model.json', 'r')
loaded_model_json = json_file.read()
#json_file.close()
loaded_model = model_from_json(loaded_model_json)
# load weights into new model
loaded_model.load_weights("best_model.hdf5")
loaded_model.compile(optimizer='adam',loss='binary_crossentropy')

from sklearn.metrics import confusion_matrix
from sklearn.metrics import plot_confusion_matrix
import seaborn as sn
#from scipy.stats import threshold

#df1 = pd.read_hdf('/content/drive/MyDrive/dataset_2500000.h5')
#y_true=np.array(df1['label'])
#features=df.iloc[:,0:24]
# Scale the train set
#scaler = StandardScaler().fit(features)
#features = scaler.transform(features)
y_pred=model.predict(X_test)
y_pr=[]
for a in y_pred:
  if a > 0.5:
    y_pr.append(1)
  else:
    y_pr.append(0)
#print(y_pr)
y_test['predicted_label']=y_pr
print(y_test)

import math
fig = plt.figure(figsize=(10,4))
ax1=fig.add_subplot(131)
#fig.suptitle('TEST MCE distribution', fontsize=12)
one_true = y_test[y_test['label']==1]['MCE']
zero_true = y_test[y_test['label']==0]['MCE']

#plt.hist(one_true,bins=100,stacked=True,label='MCE test sample, bkg')
ax1.hist([zero_true,one_true],bins=100,stacked=True,label=['signal','bkg'])
xposition = [3,4,4+math.log10(3),5]
for xc in xposition:
    ax1.axvline(x=xc, color='r', linestyle='--')
ax1.set_xlabel('log10(MCE/GeV)')
ax1.legend()
ax9=fig.add_subplot(132)
ax9.set_xlabel('log10(MCE/GeV)')
ax9.hist(one_true,bins=100,alpha=0.7)
ax9.set_title('(a)nu_mu_CC_MCE_Distro')
ax10=fig.add_subplot(133)
ax10.set_xlabel('log10(MCE/GeV)')
ax10.hist(zero_true,bins=100,alpha=0.7)
ax10.set_title('(a)nu_l_CC/NC_MCE_Distro')
plt.tight_layout()

df_till_1TeV=y_test[y_test['MCE']<3]
y_te=df_till_1TeV['label']
y_pr=df_till_1TeV['predicted_label']
df_cm=confusion_matrix(y_te,y_pr,normalize='true')
print(df_cm)
ax = plt.axes()
ax.set_title('Below 1 TeV')
sn.heatmap(df_cm, annot=True,cmap=plt.cm.Blues, fmt='g',ax=ax)
ax.set_ylabel('true')
ax.set_xlabel('predicted')

df_1TeV_10TeV=y_test[y_test['MCE'].between(3,4)]
y_te1=df_1TeV_10TeV['label']
y_pr1=df_1TeV_10TeV['predicted_label']
df_cm1=confusion_matrix(y_te1,y_pr1,normalize='true')
ax1 = plt.axes()
ax1.set_title('1 - 10 TeV')
sn.heatmap(df_cm1, annot=True,cmap=plt.cm.Blues, fmt='g',ax=ax1)
ax1.set_ylabel('true')
ax1.set_xlabel('predicted')

import math
df_10TeV_30TeV=y_test[y_test['MCE'].between(4,4+math.log10(3))]
y_te2=df_10TeV_30TeV['label']
y_pr2=df_10TeV_30TeV['predicted_label']
df_cm2=confusion_matrix(y_te2,y_pr2,normalize='true')
ax2 = plt.axes()
ax2.set_title('10 - 30 TeV')
sn.heatmap(df_cm2, annot=True,cmap=plt.cm.Blues, fmt='g',ax=ax2)
ax2.set_ylabel('true')
ax2.set_xlabel('predicted')

import math
df_30TeV_100TeV=y_test[y_test['MCE'].between(4+math.log10(3),5)]
y_te3=df_30TeV_100TeV['label']
y_pr3=df_30TeV_100TeV['predicted_label']
df_cm3=confusion_matrix(y_te3,y_pr3,normalize='true')
ax3 = plt.axes()
ax3.set_title('30 - 100 TeV')
sn.heatmap(df_cm3, annot=True,cmap=plt.cm.Blues, fmt='g',ax=ax3)
ax3.set_ylabel('true')
ax3.set_xlabel('predicted')

import math
df_above_100TeV=y_test[y_test['MCE']>5]
y_te4=df_above_100TeV['label']
y_pr4=df_above_100TeV['predicted_label']
df_cm4=confusion_matrix(y_te4,y_pr4,normalize='true')
ax4 = plt.axes()
ax4.set_title('Above 100 TeV')
sn.heatmap(df_cm4, annot=True,cmap=plt.cm.Blues, fmt='g',ax=ax4)
ax4.set_ylabel('true')
ax4.set_xlabel('predicted')

from sklearn.metrics import roc_curve
rcyt,rcyp, _ = roc_curve(y_te,y_pr)
rcyt1,rcyp1, _ = roc_curve(y_te1,y_pr1)
rcyt2,rcyp2, _ = roc_curve(y_te2,y_pr2)
rcyt3,rcyp3, _ = roc_curve(y_te3,y_pr3)
rcyt4,rcyp4, _ = roc_curve(y_te4,y_pr4)

ax5 = plt.axes()

ax5.plot(rcyt, rcyp, label='Below 1 TeV')
ax5.plot(rcyt1, rcyp1, label='1 -10 TeV')
ax5.plot(rcyt2, rcyp2, label='10 - 30 TeV')
ax5.plot(rcyt3, rcyp3, label='30 - 100 TeV')
ax5.plot(rcyt4, rcyp4, label='above 100 TeV')
ax5.legend()
plt.title('ROC_CURVE')

tp=np.array([df_cm[0][0],df_cm1[0][0],df_cm2[0][0],df_cm3[0][0],df_cm4[0][0]])*100
fp=np.array([df_cm[0][1],df_cm1[0][1],df_cm2[0][1],df_cm3[0][1],df_cm4[0][1]])*100
fn=np.array([df_cm[1][0],df_cm1[1][0],df_cm2[1][0],df_cm3[1][0],df_cm4[1][0]])*100
tn=np.array([df_cm[1][1],df_cm1[1][1],df_cm2[1][1],df_cm3[1][1],df_cm4[1][1]])*100
energy=[1.5,3.5,4+math.log10(2),3+math.log10(65),7]
ax6 = plt.axes()
ax6.plot(energy,tp,label='TP (0-0)')
ax6.plot(energy,fp,label='FP (0-1)')
ax6.plot(energy,fn,label='FN (1-0)')
ax6.plot(energy,tn,label='TN (1-1)')
ax6.plot(energy,tp,'ok')
ax6.plot(energy,fp,'ok')
ax6.plot(energy,fn,'ok')
ax6.plot(energy,tn,'ok')
ax6.legend()
ax6.set_ylabel('percentage(%)')
ax6.set_xlabel('log_10(MCE/GeV)')