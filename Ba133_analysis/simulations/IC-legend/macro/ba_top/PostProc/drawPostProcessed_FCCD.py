# draw energy histograms from each detector

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

    hdf5_path = "/lfs/l1/legend/users/aalexander/hdf5_output/processed/"

    if(len(sys.argv) != 2):
        print('Usage: drawPostProcessed_FCCD.py [e.g. IC160A_ba_top_coll_01]')
        sys.exit()

    MC_file_id = sys.argv[1] #MC_file_id = 'IC160A_ba_top_coll_01'
    MC_file = hdf5_path+"processed_detector_"+MC_file_id+'.hdf5'

    binwidth = 0.15 #keV
    #____________PLOT and Fit _________

    #no FCCD
    # plt.figure()
    # counts, bins = get_histo_energies(MC_file, binwidth)
    # plt.xlabel("Energy [keV]")
    # plt.ylabel("Counts")
    # plt.xlim(0, 450)
    # plt.yscale("log")
    # plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/simulations/IC-legend/macro/ba_top/PostProc/plots/"+MC_file_id+'.png')

    
    xmin_356, xmax_356 = 350, 360.5 #362 #360.5 for gammas #2 #360 #kev 
    plt.figure()
    counts, bins = get_histo_energies(MC_file, binwidth)
    popt, pcov, xfit = fit_peak_356_2("Energy (keV)", bins, counts, xmin_356, xmax_356)
    a,b,c,d,e = popt[0],popt[1],popt[2],popt[3],popt[4]
    plt.xlim(xmin_356, xmax_356) 
    plt.ylim(10, 10**7)
    plt.yscale("log")
    plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/simulations/IC-legend/macro/ba_top/PostProc/plots/"+MC_file_id+'_356keV.png')

    C_356, C_356_err = gauss_count(a, c, np.sqrt(pcov[0][0]), np.sqrt(pcov[2][2]), binwidth)
    print("gauss count 356keV: ", C_356, " +/- ", C_356_err )


    xmin_81, xmax_81 = 77, 84
    plt.figure()
    counts, bins = get_histo_energies(MC_file, binwidth)
    popt, pcov, xfit = fit_double_peak_81("Energy (keV)", bins, counts, xmin_81, xmax_81)
    a,b,c,d,e,f,g,h = popt[0],popt[1],popt[2],popt[3],popt[4],popt[5],popt[6],popt[7] 
    plt.xlim(xmin_81, xmax_81) 
    plt.ylim(5*10**2, 5*10**6)
    plt.yscale("log")
    plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/simulations/IC-legend/macro/ba_top/PostProc/plots/"+MC_file_id+'_81keV.png')

    R = 2.65/32.9
    C_81, C_81_err = gauss_count(a, c, np.sqrt(pcov[0][0]), np.sqrt(pcov[2][2]), binwidth)
    C_79, C_79_err = gauss_count(R*a, e, R*np.sqrt(pcov[0][0]), np.sqrt(pcov[4][4]), binwidth)
    print("gauss count 81: ", C_81, " +/- ", C_81_err )
    print("gauss count 79.6: ", C_79, " +/- ", C_79_err )

    print("")
    O_Ba133 = (C_79 + C_81)/C_356
    O_Ba133_err = O_Ba133*np.sqrt((C_79_err**2 + C_81_err**2)/(C_79+C_81)**2 + (C_356_err/C_356)**2)
    print("O_BA133 = " , O_Ba133, " +/- ", O_Ba133_err)


    #Save count values to json file
    dlt_observables = {
        "C_356": C_356,
        "C_356_err" : C_356_err,
        "C_81" : C_81,
        "C_81_err" : C_81_err,
        "C_79" : C_79,
        "C_79_err" : C_79_err,
        "O_Ba133" : O_Ba133,
        "O_Ba133_err" : O_Ba133_err
    }
    with open("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/simulations/IC-legend/macro/ba_top/PostProc/"+MC_file_id+"_dlt_observables.json", "w") as outfile: 
        json.dump(dlt_observables, outfile)


    # #_____________PROCESS AND PLOT FCCDS_____________

    # # #graph of all FCCDs
    # # FCCD_list = ['none', 0.25, 0.5, 0.75, 1.0, 1.25, 1.5] #mm
    # # process_FCCDs(MC_file_id, FCCD_list, binwidth)
    # # plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/simulations/IC-legend/macro/ba_top/PostProc/plots/FCCDs_"+MC_file_id+'.png')
    
    # # #repeat for just none and 1.5 - FCCDs2
    # # FCCD_list = ['none', 1.5] #mm, FCCDS2
    # # process_FCCDs(MC_file_id, FCCD_list, binwidth)
    # # plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/simulations/IC-legend/macro/ba_top/PostProc/plots/FCCDs2_"+MC_file_id+'.png')


    #plt.show() 
    
    print("done")


def get_histo_energies(MC_file, binwidth):

    df =  pd.read_hdf(MC_file, key="procdf")
    energies = df['energy']
    energies = energies*1000

    #binwidth = 0.15
    
    counts, bins, bars = plt.hist(energies, bins = np.arange(min(energies), max(energies) + binwidth, binwidth), histtype = 'step') #, linewidth = '0.35')
    no_events = energies.size #=sum(counts)
    print("No. events: ", no_events) 

    return counts, bins


def process_FCCDs(MC_file_id, FCCD_list, binwidth):
    "process and plot for a list of different FCCDs"

    hdf5_path = "/lfs/l1/legend/users/aalexander/hdf5_output/processed/"

    fig, ax = plt.subplots()
    #binwidth = 0.15 #5 #0.1 kev = rough min resolution
    
    for FCCD in FCCD_list:

        if FCCD == 'none':
            df =  pd.read_hdf(hdf5_path+"processed_detector_"+MC_file_id+'.hdf5', key="procdf")
        else: 
            df =  pd.read_hdf(hdf5_path+"processed_FCCD"+str(FCCD)+"mm_detector_"+MC_file_id+'.hdf5', key="procdf")
        
        energies = df['energy']
        energies = energies*1000
        counts, bins, bars = plt.hist(energies, bins = np.arange(min(energies), max(energies) + binwidth, binwidth), label = 'FCCD: '+str(FCCD)+' mm', histtype = 'step', linewidth = '0.35')
        no_events = energies.size #=sum(counts)
        print("No. events: ", no_events) 

    plt.xlabel("Energy [keV]")
    plt.ylabel("Counts")
    plt.xlim(0, 450)
    plt.yscale("log")
    #plt.ylim(1,10**7)  
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    info_str = r'binwidth = $%.2f$ keV' % (binwidth)
    ax.text(0.67, 0.97, info_str, transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props)
    plt.legend(loc='lower left', fontsize = 7.5)
    #plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/simulations/IC-legend/macro/ba_top/PostProc/plots/FCCDs_"+MC_file_id+'.png')
    




if __name__ == "__main__":
    main()