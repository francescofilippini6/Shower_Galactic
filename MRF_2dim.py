import ROOT
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
ROOT.gInterpreter.ProcessLine('#include <TFeldmanCousins.h>')
ROOT.gInterpreter.Declare("""
Double_t MeanUL(Double_t Nbackground)
{
  Double_t CL=0.9; // set confidence level
  TFeldmanCousins FC(CL);
  FC.SetMuMax(100.0); // maximum value of signal to calculate the tables!
                      // increase it for greater values of Nbackground!
  FC.SetMuStep(0.1);
  Double_t PoisProb=0., mu90=0., ul=0.;//, ll=0.;

  //for (Int_t Nobserved=0; Nobserved<200; Nobserved++) // decrease this value for Nbackground<20
  //for (Int_t Nobserved=0; Nobserved<100; Nobserved++) // decrease this value for Nbackground<20
  /* * * * * * * * * * NON FUNZIONA CORRETTAMENTE PER N>20 * * * * * * * * * */
 for (Int_t Nobserved=0; Nobserved<200; Nobserved++) // increase this value for Nbackground>20
    {
      ul = FC.CalculateUpperLimit(Nobserved, Nbackground);
      //ll = FC.GetLowerLimit();
      PoisProb = TMath::Poisson(Nobserved,Nbackground);

      mu90 = mu90 + (ul * PoisProb);
      //      cout << "\t" << Nobserved << "\t" << mu90 << endl; // DEBUG!
    }

  return mu90;
}
""")

#def mrfplotter(df):
#    fig=plt.figure()
#    ax=fig.add_subplot(111)
#    x=df['center_bin']
#    y= np.linspace(0.1,0.9,9, endpoint=True)
#    X,Y = plt.meshgrid(x, y) # grid of point
#    Z = (X, Y) # evaluation of the function on the grid
#    
#    im = imshow(Z,cmap=cm.RdBu) # drawing the function
#    # adding the Contour lines with labels
#    cset = contour(Z,arange(-1,1.5,0.2),linewidths=2,cmap=cm.Set2)
#    clabel(cset,inline=True,fmt='%1.1f',fontsize=10)
#    colorbar(im) # adding the colobar on the right
#    ax.plot(df['center_bin'],df['mrf'],'+g')
#    ax.set_yscale('log')
#    ax.set_ylabel('MRF')
#    ax.set_xlabel(r'$\log_{10}(E_{Tantra}/(GeV))$')
#    plt.show()
#    return fig


def mrf_eneergy_dimension(df,lower_index,upper_index):
    MRF=[]
    for a in range(lower_index,upper_index):
        if df['OffEntries'][a]==0:
            MRF.append(0)
        elif df['OnEntries'][a]==0:
            MRF.append(0)
        else:
            print(df['OffEntries'][a])
            fcul=ROOT.MeanUL(df['OffEntries'][a])
            ns=df['OnEntries'][a]
            MRF.append(fcul/ns)
    #df['mrf']=MRF
    #return df
    return MRF
    
if __name__ == "__main__":
    #df=pd.read_csv('ON_OFF_histo.csv')
    df=pd.read_csv('MRF_2dim.csv')
    print(df)
    #listofdf=[]
    #----------------------------------------------
    # initializing the 2d istogram with random values
    #----------------------------------------------
    xedges = np.linspace(0,0.9,10, endpoint=True)
    dfhh=df[0:50]
    energy_edge=np.array(dfhh['bin_edges'])
    #print(energy_edge)
    #print(len(energy_edge))
    yedges=np.insert(energy_edge,0, -0.041204)
    #print(yedges)
    #print(len(yedges))
    x = np.random.normal(0.5, 0.3, 100)
    y = np.random.normal(1, 3, 100)
    H, xedges, yedges = np.histogram2d(x,y, bins=(xedges, yedges))
    #print(len(H[8]))
    aa=[0,50,100,150,200,250,300,350,400,450]
    for i in range(9):
        print(aa[i],aa[i+1])
        df1=df[aa[i]:aa[i+1]]
        H[i]=mrf_eneergy_dimension(df1,aa[i],aa[i+1])
        print("histogram",H[i])
        
    myextent  =[xedges[0],xedges[-1],yedges[0],yedges[-1]]
    
    fig=plt.figure()
    ax=fig.add_subplot(121)
    AA=ax.imshow(H.T,origin='lower',extent=myextent,interpolation='nearest',aspect='auto')
    ax.set_title('MRF')
    ax.set_ylabel(r'$\log_{10}(E_{Tantra}/(GeV))$')
    ax.set_xlabel('ANN prediciton')
    fig.colorbar(AA, ax=ax)
    ax1=fig.add_subplot(122)
    aaa=ax1.contour(H.T,extent=myextent,linewidths=1,cmap='Set2')
    #ax1.clabel(aaa,inline=True,fmt='%1.1f',fontsize=6)
    plt.show()
    
    my_df = pd.DataFrame(H)
    my_df.to_csv('my_MRF_2dim.csv', index=False, header=False)
    
    
