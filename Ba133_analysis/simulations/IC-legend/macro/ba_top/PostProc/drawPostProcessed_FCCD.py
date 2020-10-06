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

#CURRENT CODE - plot spectra, fit, compare with data and deadlayer post processing (currently testing)

def main(): 


    #print date and time for log:
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S") # dd/mm/YY H:M:S
    print("")
    print("date and time =", dt_string)	
    print("")

    hdf5_path = "/lfs/l1/legend/users/aalexander/hdf5_output/processed/"

    if(len(sys.argv) != 2):
        print('Usage: drawPostProcessed_FCCD.py [e.g. IC160A_ba_top_81mmNEW4_01]')
        sys.exit()

    MC_file_id = sys.argv[1] #MC_file_id = 'IC160A_ba_top_coll_01'
    MC_file = hdf5_path+"processed_detector_"+MC_file_id+'.hdf5'

    binwidth = 0.15 #keV

    #____________PLOT Spectra _________

    print("plotting whole simulated spectrum...")
    
    # plot full spectrum

    df =  pd.read_hdf(MC_file, key="procdf")
    energies = df['energy']
    energies = energies*1000
    no_events = energies.size #=sum(counts)
    print("No. events: ", no_events) 
    bins = np.arange(min(energies), max(energies) + binwidth, binwidth)

    fig, ax = plt.subplots()
    counts, bins, bars = plt.hist(energies, bins = bins) #, linewidth = '0.35')
    plt.xlabel("Energy [keV]")
    plt.ylabel("Counts")
    plt.xlim(0, 450)
    plt.yscale("log")
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    info_str = '\n'.join((r'# events = $%.0f$' % (no_events), r'binwidth = $%.2f$ keV' % (binwidth)))
    ax.text(0.67, 0.97, info_str, transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props)
    plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/simulations/IC-legend/macro/ba_top/PostProc/plots/"+MC_file_id+'.png')


    # #________Fit peaks of interest_______
    print("")
    print("Fitting peaks of interest...")

    xmin_356, xmax_356 = 350, 362 #362 #360.5 for gammas #2 #360 #kev 
    plt.figure()
    counts, bins, bars = plt.hist(energies, bins = bins, histtype = 'step') #, linewidth = '0.35')
    popt, pcov, xfit = fit_peak_356_2("Energy (keV)", bins, counts, xmin_356, xmax_356)
    a,b,c,d,e = popt[0],popt[1],popt[2],popt[3],popt[4]
    amplitude356_sim = gaussian_and_bkg_2(b, a, b, c, d, e)
    plt.xlim(xmin_356, xmax_356) 
    plt.ylim(10, 10**7)
    plt.yscale("log")
    plt.xlabel("Energy [keV]")
    plt.ylabel("Counts")
    plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/simulations/IC-legend/macro/ba_top/PostProc/plots/"+MC_file_id+'_356keV.png')

    C_356, C_356_err = gauss_count(a, c, np.sqrt(pcov[0][0]), np.sqrt(pcov[2][2]), binwidth)
    print("gauss count 356keV: ", C_356, " +/- ", C_356_err )


    xmin_81, xmax_81 = 77, 84
    plt.figure()
    counts, bins, bars = plt.hist(energies, bins = bins, histtype = 'step') #, linewidth = '0.35')
    popt, pcov, xfit = fit_double_peak_81("Energy (keV)", bins, counts, xmin_81, xmax_81)
    a,b,c,d,e,f,g,h = popt[0],popt[1],popt[2],popt[3],popt[4],popt[5],popt[6],popt[7] 
    plt.xlim(xmin_81, xmax_81) 
    plt.ylim(5*10**2, 5*10**6)
    #plt.ylim(10**3, 10**7) #gammas_81mmNEW
    plt.yscale("log")
    plt.xlabel("Energy [keV]")
    plt.ylabel("Counts")
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


    # #__________compare against real data__________
    print("")
    print("plotting simulation against actual data...")

    #this code below is from Ba133_dlt_analysis.py
    detector = "I02160A" #read tier 2 runs for Ba data
    t2_folder = "/lfs/l1/legend/detector_char/enr/hades/char_data/"+detector+"/tier2/ba_HS4_top_dlt/pygama/"
    keys, data = read_all_t2(t2_folder)
    data_size = data.size #all events
    key_data = obtain_key_data(data, keys, "e_ftp", data_size)
    with open('/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/data/calibration_coef.json') as json_file:
        calibration_coefs = json.load(json_file)
        m = calibration_coefs['m']
        m_err = calibration_coefs['m_err']
        c = calibration_coefs['c']
        c_err = calibration_coefs['c_err']
        # a_quad = calibration_coefs['a_quad']
        # a_quad_err = calibration_coefs['a_quad_err']
        # b_quad = calibration_coefs['b_quad']
        # b_quad_err = calibration_coefs['b_quad_err']
        # c_quad = calibration_coefs['c_quad']
        # c_quad_err = calibration_coefs['c_quad_err']

    calibrated_energy_data = (key_data-c)/m

    #plot absolutes
    bins_data = bins = np.arange(0, 450, binwidth)
    fig, ax = plt.subplots()
    counts_data, bins_cal_data, bars_data = plt.hist(calibrated_energy_data, bins=bins_data, label = "data", histtype = 'step', linewidth = '0.35')
    counts, bins, bars = plt.hist(energies, bins = bins, label = "MC: No FCCD", histtype = 'step', linewidth = '0.35')
    #counts, bins, bars = plt.hist(energies, bins = bins, label = "MC: FCCD 0.568mm", histtype = 'step', linewidth = '0.35')
    plt.xlabel("Energy [keV]")
    plt.ylabel("Counts")
    plt.xlim(0, 450)
    plt.yscale("log")
    #props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    #info_str = '\n'.join((r'# events = $%.0f$' % (no_events), r'binwidth = $%.2f$ keV' % (binwidth)))
    #ax.text(0.67, 0.97, info_str, transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props)
    plt.legend(loc = "lower left")
    plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/simulations/IC-legend/macro/ba_top/PostProc/plots/"+MC_file_id+'_DATAcomparison.png')

    #scale up data to same amplitude 356 peak as simulation
    amplitude356_data = gaussian_and_bkg_2(356, 4.6*10**4, 356, 0.423, 2.64, 205) #from old plots
    print("amplitude 356keV data: ", amplitude356_data)
    print("amplitude 356keV simulation: ", amplitude356_sim)
    R_simdata_356 = amplitude356_sim/amplitude356_data #ratio of sim to data for 356 peak

    fig, ax = plt.subplots()
    counts_data, bins_cal_data, bars_data = plt.hist(calibrated_energy_data, bins=bins_data, weights=R_simdata_356*np.ones_like(calibrated_energy_data),  label = "data (scaled)", histtype = 'step', linewidth = '0.35')
    counts, bins, bars = plt.hist(energies, bins = bins, label = "MC: No FCCD", histtype = 'step', linewidth = '0.35')
    #counts, bins, bars = plt.hist(energies, bins = bins, label = "MC:FCCD 0.568mm", histtype = 'step', linewidth = '0.35')
    plt.xlabel("Energy [keV]")
    plt.ylabel("Counts")
    plt.xlim(0, 450)
    plt.yscale("log")
    #props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    #info_str = '\n'.join((r'# events = $%.0f$' % (no_events), r'binwidth = $%.2f$ keV' % (binwidth)))
    #ax.text(0.67, 0.97, info_str, transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props)
    plt.legend(loc = "lower left")
    plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/simulations/IC-legend/macro/ba_top/PostProc/plots/"+MC_file_id+'_DATAcomparison_scaled.png')

    #calculate DATA/MC for each energy bin and export
    print("")
    print("calculating data/MC ratios...")

    Data_MC_ratios = []
    print(len(counts_data))
    print(len(bins_data))
    for index, bin in enumerate(bins_data[1:]):
        data = counts_data[index]
        MC = counts[index]
        if MC == 0:
            ratio = 0.
        else:
            try: 
                ratio = data/MC
            except:
                ratio = 0 #if MC=0 and dividing by 0
        Data_MC_ratios.append(ratio)

    Data_MC_ratios_df = pd.DataFrame({'ratio': Data_MC_ratios})
    Data_MC_ratios_df['bin'] = bins_data[1:]
    print(Data_MC_ratios_df)
    Data_MC_ratios_df.to_csv("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/simulations/IC-legend/macro/ba_top/PostProc/"+MC_file_id+"_DataMCRatios.csv", index=False)
    




    # #_____________PROCESS AND PLOT FCCDS_____________ #in development/testing

    # # #graph of all FCCDs
    # # FCCD_list = ['none', 0.25, 0.5, 0.75, 1.0, 1.25, 1.5] #mm
    # # process_FCCDs(MC_file_id, FCCD_list, binwidth)
    # # plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/simulations/IC-legend/macro/ba_top/PostProc/plots/FCCDs_"+MC_file_id+'.png')
    
    # # #repeat for just none and 1.5 - FCCDs2
    # # FCCD_list = ['none', 1.5] #mm, FCCDS2
    # # process_FCCDs(MC_file_id, FCCD_list, binwidth)
    # # plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/simulations/IC-legend/macro/ba_top/PostProc/plots/FCCDs2_"+MC_file_id+'.png') 
    
    print("done")


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