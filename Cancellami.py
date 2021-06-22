import sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns
import numpy as np

    
if __name__ == "__main__":
    #df=pd.read_csv('ON_OFF_histo.csv')
    df=pd.read_hdf('MRF_3D.h5')
    print(df.iloc[1,1])
    df.to_csv('MRF_3D.csv')
    #    listofdf=[]
#    #----------------------------------------------
#    # initializing the 2d istogram with random values
#    #----------------------------------------------
#    yedges = np.linspace(0.025,0.925,19, endpoint=True)
#    xedges= np.linspace(-0.075,0.575,14, endpoint=True)#-0.35,12
#    x = np.random.normal(0.5, 0.3, 100)
#    y = np.random.normal(1, 3, 100)
#    H, xedges, yedges = np.histogram2d(x,y, bins=(xedges, yedges))
#    print(H[0])
#    counter=0
#    for i in range(13):
#        df1=df[counter:counter+3]
#        print(df1)
#        print('\n')
#        off=np.array(df1.iloc[0])
#        on=np.array(df1.iloc[2])
#        counter+=3
#        MRF=[]
#        minima=[]
#        for a in range(len(on)):
#            MRF.append(mrf_2dimension(off[a],on[a]))
#        #print(on)
#        minima.append(np.argmin(MRF))
#        H[i] = MRF
#
#    print(minima)
#    #print(H)
#    myextent  =[xedges[0],xedges[-1],yedges[0],yedges[-1]]
#    fig=plt.figure()
#    ax=fig.add_subplot(121)
#    AA=ax.imshow(H.T,origin='lower',extent=myextent,interpolation='nearest',aspect='auto',norm = LogNorm())
#    ax.set_title('MRF')
#    ax.set_ylabel('NN output')
#    ax.set_xlabel('BDT cut')
#    fig.colorbar(AA, ax=ax)
#    ax1=fig.add_subplot(122)
#    aaa=ax1.contourf(H.T,extent=myextent,linewidths=1,cmap='Set2',norm = LogNorm())
#    #ax1.clabel(aaa,inline=True,fmt='%1.1f',fontsize=6)
#    plt.show()
#    #my_df = pd.DataFrame(H)
#    #my_df.to_csv('my_MRF_2dim.csv', index=False, header=False)
#    
#    
#
