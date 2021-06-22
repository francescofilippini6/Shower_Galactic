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
  FC.SetMuMax(100.0); //500 maximum value of signal to calculate the tables!
                      // increase it for greater values of Nbackground!
  FC.SetMuStep(0.1); //0.1 the step in signal to use when generating tables
  Double_t PoisProb=0., mu90=0., ul=0.;//, ll=0.;

  //for (Int_t Nobserved=0; Nobserved<200; Nobserved++) // decrease this value for Nbackground<20
  //for (Int_t Nobserved=0; Nobserved<100; Nobserved++) // decrease this value for Nbackground<20
  /* * * * * * * * * * NON FUNZIONA CORRETTAMENTE PER N>20 * * * * * * * * * */
 for (Int_t Nobserved=0; Nobserved<200; Nobserved++) //1000 before increase this value for Nbackground>20
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


def mrf_3dimension(bkg,sign):
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
    df=pd.read_csv('MRF_3D.csv')
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    print(df)
    bdt_bin=np.linspace(-0.05,0.55,13,endpoint=True)   # step 0.05
    ann_bin=np.linspace(0.05,0.9,18, endpoint=True)
    energy_bin=np.linspace(0.76,8.08,46,endpoint=True) # max value 9.06
    table = pd.DataFrame(index = ['-0.05','0.0','0.05','0.10','0.15','0.20','0.25','0.3' ,'0.35','0.4','0.45','0.5','0.55'],columns = ['0.05','0.10','0.15','0.20','0.25','0.3' ,'0.35','0.4','0.45','0.5','0.55','0.6' ,'0.65','0.7', '0.75','0.8' ,'0.85','0.9'])
    for row,bdt_cut in enumerate(bdt_bin):
        for column,ann_cut in enumerate(ann_bin):
            print("ann cut:",ann_cut)
            aa=df.iloc[row,column]
            #print(aa)
            aa=aa.replace("[", "")
            aa=aa.replace("]", "")
            ee=aa.split(',')
            ae=[]
            for a in ee:
                ae.append(float(a))
            off=ae[:50]
            on=ae[51:100]
            cosmic=ae[101:150]
            MRF=[]
            for a in range(len(on)):
                MRF.append(mrf_3dimension(off[a],cosmic[a]))
            table.iloc[row,column]=MRF
            print(MRF)
    print("saving file")
    table.to_csv('MRF_cut_energy.csv')
    
            
            
