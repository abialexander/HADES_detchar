# draw energy histograms from each detector

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#plt.style.use("mplstyle.txt")


def main(): 

    hdf5_path = "/lfs/l1/legend/users/aalexander/hdf5_output/processed/"

    # if(len(sys.argv) != 2):
    #     print('Usage: drawPostProcessedHdf5.py [detector_IC160A_ba_top_coll_01]')
    #     sys.exit()

    # MC_file = sys.argv[1]

    MC_file = 'detector_IC160A_ba_top_coll_01'

    FCCD_list = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5] #mm
    #FCCD_list = [ 1.5] #mm, FCCDS2

    fig, ax = plt.subplots()
    binwidth = 0.15 #5 #0.1 kev = rough min resolution

    FCCD = "none"
    df =  pd.read_hdf(hdf5_path+"processed_"+MC_file+'.hdf5', key="procdf")
    energies = df['energy']
    energies = energies*1000
    counts, bins, bars = plt.hist(energies, bins = np.arange(min(energies), max(energies) + binwidth, binwidth), label = 'FCCD: '+FCCD, histtype = 'step')
    no_events = energies.size #=sum(counts)
    print("No. events: ", no_events)

    for FCCD in FCCD_list:
        df =  pd.read_hdf(hdf5_path+"processed_FCCD"+str(FCCD)+"mm_"+MC_file+'.hdf5', key="procdf")
        energies = df['energy']
        energies = energies*1000
        counts, bins, bars = plt.hist(energies, bins = np.arange(min(energies), max(energies) + binwidth, binwidth), label = 'FCCD: '+str(FCCD)+' mm', histtype = 'step')
        no_events = energies.size #=sum(counts)
        print("No. events: ", no_events) 

    plt.xlabel("Energy [keV]")
    plt.ylabel("Counts")
    plt.xlim(0, 450)
    plt.yscale("log")
    plt.ylim(1,10**4)  
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    info_str = r'binwidth = $%.2f$ keV' % (binwidth)
    #info_str = '\n'.join((r'# events = $%.0f$' % (no_events), r'binwidth = $%.2f$ keV' % (binwidth)))
    ax.text(0.67, 0.97, info_str, transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props)
    plt.legend(loc='upper left')

    plt.savefig("plots/FCCDs_"+MC_file+'.png')



if __name__ == "__main__":
    main()