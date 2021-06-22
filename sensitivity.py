import sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns
import numpy as np
    
if __name__ == "__main__":
    #df=pd.read_csv('ON_OFF_histo.csv')
    df=pd.read_csv('Calc_MRF_3D_step_smaller.csv')
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    #print(df.head(10))
    #print(df['0.5'][1])
    bdt_bin=np.linspace(-0.05,0.55,13,endpoint=True)   # step 0.05
    ann_bin=np.linspace(0.05,0.9,18, endpoint=True)
    energy_bin=np.linspace(0.76,9.06,51,endpoint=True)
    minima=[]
    for row,bdt_cut in enumerate(bdt_bin):
        for column,ann_cut in enumerate(ann_bin):
            aa=df.iloc[row,column]
            #print(aa)
            aa=aa.replace("[", "")
            aa=aa.replace("]", "")
            ee=aa.split(',')
            #ae=[]
            for ene,a in enumerate(ee):
                if float(a) == 0.0:
                    minima.append(100000)
                else:
                    minima.append(float(a))
    final=min(minima)
    print("MRF final value:", final)
    
    for row,bdt_cut in enumerate(bdt_bin):
        for column,ann_cut in enumerate(ann_bin):
            aab=df.iloc[row,column]
            aab=aab.replace("[", "")
            aab=aab.replace("]", "")
            eeb=aab.split(',')
            aeb=[]
            for ener,ab in enumerate(eeb):
                if float(ab) == final:
                    print('ANN, BDT, ENE: ',ann_cut,bdt_cut,energy_bin[ene])
                else:
                    continue
