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
    signal=[4.1705,4.0088,3.57412651430477,3.40530226216786,3.011526845293326,2.875068785046655,2.295507792360021]
    bkg=[9.50,7.83,5.347275129578436,5.0624569411158635,3.47393009675826,3.295897313660634,2.1119067362014206]
    for i in range(len(signal)):
        off=bkg[i]
        on=signal[i]
        MRF=mrf_2dimension(off,on)
        print(MRF)
        
    
    
    
    
