#plot data/MC ratios for different FCCDs
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#plt.style.use("mplstyle.txt")
from datetime import datetime

#import fitting functions
import sys
sys.path.append('/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/data/')
from Ba133_dlt_analysis import * 

def main():

    #print date and time for log:
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S") # dd/mm/YY H:M:S
    print("")
    print("date and time =", dt_string)	
    print("")


    ratios_noFCCD_df = pd.read_csv("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/simulations/IC-legend/macro/ba_top/PostProc/IC160A_ba_top_81mmNEW4_01_DataMCRatios.csv")
    ratios_noFCCD = ratios_noFCCD_df['ratio']
    energies = ratios_noFCCD_df['bin']

    ratios_FCCD0568mm_df = pd.read_csv("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/simulations/IC-legend/macro/ba_top/PostProc/IC160A_ba_top_81mmNEW4_01_FCCD0.568mm_DataMCRatios.csv")
    ratios_FCCD0568mm = ratios_FCCD0568mm_df['ratio']
    
    
    #ratios_FCCD0568mm.to_numpy()

    plt.figure()
    plt.plot(energies, ratios_noFCCD, 'o', ms=1.25,label = "no FCCD")
    plt.plot(energies, ratios_FCCD0568mm, 'o', ms =1.25, label = "FCCD: 0.568mm")
    ones = [1]*len(energies)
    plt.plot(energies, ones, "k-.")
    plt.xlabel('Energy (keV)')
    plt.ylabel('Data/MC')
    plt.yscale("log")
    plt.ylim(0,50)
    plt.xlim(0,450)
    plt.legend()
    plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/simulations/IC-legend/macro/ba_top/PostProc/plots/Data_MC_ratios.png")




if __name__ == "__main__":
    main()