# draw energy histograms from each detector

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#plt.style.use("mplstyle.txt")


def main(): 

    hdf5_path = "/lfs/l1/legend/users/aalexander/hdf5_output/processed/"

    if(len(sys.argv) != 2):
        print('Usage: drawPostProcessedHdf5.py [processed.hdf5]')
        sys.exit()

    filename = sys.argv[1]
    plotname = filename.strip('.hdf5')

    df =  pd.read_hdf(hdf5_path+filename, key="procdf")

    # bins = np.arange(0, 4000, 1) #want rougly binwidth to be det resolution, e.g. 0.1keV 
    # plt.figure()
    # for det, detdf in df.groupby('detID'):
    #     (detdf["energy"]*1000).hist(bins=bins, histtype="step", label="detector {}".format(det))
    # plt.xlabel("Energy [keV]")
    # plt.ylabel("Counts")
    # plt.gca().set_xlim(0, 450)
    # plt.gca().set_ylim(10, plt.gca().get_ylim()[1])
    # plt.gca().grid(False)
    # plt.yscale("log")
    #plt.legend(frameon=False, loc='upper right')
    #plt.savefig('example.png')
    #plt.show()


    binwidth = 0.15 #5 #0.1 kev = rough min resolution
    energies = df['energy']
    energies = energies*1000

    # #plt.figure()
    fig, ax = plt.subplots()
    #counts, bins, bars = plt.hist(energies, bins = np.arange(min(energies), max(energies) + binwidth, binwidth)) #, histtype = 'step')
    plt.hist(energies, bins = np.arange(min(energies), max(energies) + binwidth, binwidth)) #, histtype = 'step')

    

    # no_events = sum(counts)
    # print("No. events: ", no_events)

    # no_events = sum(counts)/binwidth
    # print("No. events: ", no_events)

    # no_events = sum(counts)*binwidth
    # print("No. events: ", no_events)

    # no_events = energies.size
    # print("No. events: ", no_events) #this is the correct one


    # plt.xlabel("Energy [keV]")
    # plt.ylabel("Counts")
    # plt.xlim(0, 450)
    # plt.yscale("log")
    # plt.ylim(1,10**6) 
    # props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    # info_str = r'# events = $%.0f$' % (no_events)
    # ax.text(0.67, 0.97, info_str, transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props)
    # plt.savefig("plots/"+plotname+'.png')



if __name__ == "__main__":
    main()