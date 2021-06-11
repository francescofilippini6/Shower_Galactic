import ROOT
import sys
import pandas as pd
import matplotlib.pyplot as plt
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

def mrfplotter(df):
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(df['center_bin'],df['mrf'],'+g')
    ax.set_yscale('log')
    ax.set_ylabel('MRF')
    ax.set_xlabel(r'$\log_{10}(E_{Tantra}/(GeV))$')
    plt.show()
    return fig

if __name__ == "__main__":
    #df=pd.read_csv('ON_OFF_histo.csv')
    df=pd.read_csv('Energy_Mrf.csv')
    print(df)
    MRF=[]
    for a in range(len(df['OnEntries'])):
        if df['OffEntries'][a]==0:
            MRF.append(0)
        elif df['OnEntries'][a]==0:
            MRF.append(0)
        else:
            print(df['OffEntries'][a])
            fcul=ROOT.MeanUL(df['OffEntries'][a])
            ns=df['cosmic'][a]
            MRF.append(fcul/ns)
            print(MRF)
    df['mrf']=MRF
    #    df.to_csv(r'MRF_ON_OFF_histo.csv',index=False,header=True)
    mrfplotter(df)
