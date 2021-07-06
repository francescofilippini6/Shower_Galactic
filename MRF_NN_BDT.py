import ROOT
import sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns
import numpy as np
ROOT.gInterpreter.ProcessLine('#include <TFeldmanCousins.h>')
ROOT.gInterpreter.Declare("""
Double_t MeanUL(Double_t Nbackground)
{
  Double_t CL=0.9; // set confidence level
  TFeldmanCousins FC(CL);
  FC.SetMuMax(500.0); // maximum value of signal to calculate the tables!
                      // increase it for greater values of Nbackground!
  FC.SetMuStep(1);//0.1
  Double_t PoisProb=0., mu90=0., ul=0.;//, ll=0.;

  //for (Int_t Nobserved=0; Nobserved<200; Nobserved++) // decrease this value for Nbackground<20
  //for (Int_t Nobserved=0; Nobserved<100; Nobserved++) // decrease this value for Nbackground<20
  /* * * * * * * * * * NON FUNZIONA CORRETTAMENTE PER N>20 * * * * * * * * * */
 for (Int_t Nobserved=0; Nobserved<1000; Nobserved++) //#200 before increase this value for Nbackground>20
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


def mrf_2dimension(bkg,sign):
    MRF=0
    if  bkg==0:
        MRF=0
    elif sign==0:
        MRF=0
    else:
        fcul=ROOT.MeanUL(bkg)
        MRF=fcul/sign
    return MRF
    
if __name__ == "__main__":
    #df=pd.read_csv('ON_OFF_histo.csv')
    df=pd.read_csv('cosmic_out_More_binning.csv')
    print(df)
    listofdf=[]
    #----------------------------------------------
    # initializing the 2d istogram with random values
    #----------------------------------------------
    yedges = np.linspace(0.025,0.925,19, endpoint=True)
    xedges= np.linspace(-0.075,0.575,14, endpoint=True)#-0.35,12
    x = np.random.normal(0.5, 0.3, 100)
    y = np.random.normal(1, 3, 100)
    H, xedges, yedges = np.histogram2d(x,y, bins=(xedges, yedges))
    print(H[0])
    counter=0
    #aaa=[]
    table = pd.DataFrame(index = ['-0.05','0.0','0.05','0.10','0.15','0.20','0.25','0.3' ,'0.35','0.4','0.45','0.5','0.55'],columns = ['0.05','0.10','0.15','0.20','0.25','0.3' ,'0.35','0.4','0.45','0.5','0.55','0.6' ,'0.65','0.7', '0.75','0.8' ,'0.85','0.9'])
    for i in range(13):
        df1=df[counter:counter+3]
        print(df1)
        print('\n')
        off=np.array(df1.iloc[0])
        on=np.array(df1.iloc[2])
        counter+=3
        MRF_array=[]
        #minima=[]
        for a in range(len(on)):
            MRF=mrf_2dimension(off[a],on[a])
            table.iloc[i,a]=MRF
            MRF_array.append(MRF)
        #print(on)
        #minima.append(min(MRF_array))
        H[i] = MRF
    print("saving file")
    table.to_csv('MRF_NN_BDT.csv')

    #print(minima)
    #print(H)
    myextent  =[xedges[0],xedges[-1],yedges[0],yedges[-1]]
    fig=plt.figure()
    ax=fig.add_subplot(111)
    AA=ax.imshow(H.T,origin='lower',extent=myextent,interpolation='nearest',aspect='auto',norm = LogNorm())
    ax.set_title('MRF')
    ax.set_ylabel('NN output')
    ax.set_xlabel('BDT cut')
    fig.colorbar(AA, ax=ax)
    #ax1=fig.add_subplot(122)
    #aaa=ax1.contourf(H.T,extent=myextent,linewidths=1,cmap='Set2',norm = LogNorm())
    #ax1.clabel(aaa,inline=True,fmt='%1.1f',fontsize=6)
    plt.show()
    #my_df = pd.DataFrame(H)
    #my_df.to_csv('my_MRF_2dim.csv', index=False, header=False)
    
    
