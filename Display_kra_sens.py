import sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns
import numpy as np
import csv
    
if __name__ == "__main__":
    aa=pd.read_csv('KraGamma.csv')
    bb=pd.read_csv('Fermi_point.csv')
    cc=pd.read_csv('Antares_sens.csv')
    print(aa)
    mrf=6.63527
    x=np.linspace(10**0.748,10**3,100)
    rr=np.linspace(10**-2,10**4,100)
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(aa['x'],aa['y'],label=r'$KRA_\gamma$ (cut = 5 $\times 10^7$ GeV)')
    ax.plot(x,np.array((0.145*3*mrf*4.8*10**-7)*np.array(x*1000)**-0.3),label=r'Antares sensitivity Shower ($\Gamma$=-2.3)')
    #ax.plot(x,np.array((3*0.145*2.0*10**-5)*np.array(x*1000)**-0.4),label=r'proof')
    #ax.plot(rr,np.array((3*4.8*10**-7)*np.array(rr*1000)**-0.3),label=r'Flux ($\Gamma$=-2.3)')
    ax.plot(bb['x'],bb['y'],'or',label=r'Fermi-LAT points')
    ax.plot(cc['x'],cc['y'],label=r'Antares sens 2015 ($\Gamma$=-2.4)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$E_\nu \: [TeV]$')
    ax.set_ylabel(r'$\Delta \Omega E^2_\nu \Phi^{3f} [GeV\: cm^{-2}\: s^{-1}]$')
    ax.legend()
    plt.show()
    
