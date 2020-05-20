import pandas as pd
import numpy as np
import h5py
import os
import matplotlib.pyplot as plt
import argparse
from scipy import optimize
from scipy import stats
import glob
#import pygama


def main():

    #read tier 2 runs for Ba data - I02160A
    t2_folder = "/lfs/l1/legend/detector_char/enr/hades/char_data/I02160A/tier2/ba_HS4_top_dlt/pygama/"
    keys, data = read_all_t2(t2_folder)

    no_events = data.size #all events
    key = "e_ftp"
    key_data = obtain_key_data(data, keys, key, no_events)

    plt.figure()
    counts, bins, bars = plt.hist(key_data, bins=100000)
    plt.yscale("log")
    plt.xlabel("e_ftp")
    plt.ylabel("Frequency")
    plt.xlim(0, 40000)
    plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/plots/e_ftp.png")


    #Calibrate - use linear calibration coefficients from previous work
    m, m_err = 15.786, 0.007
    c, c_err = -387.802, 11.591

    calibrated_energy = (key_data-c)/m
    plt.figure()
    counts, bins_cal, bars = plt.hist(calibrated_energy, bins=100000)
    plt.xlabel("Energy (KeV)")
    plt.ylabel("Frequency")
    plt.xlim(0, 2500)
    plt.yscale("log")
    plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/plots/calibrated_energy.png") 









def read_all_t2(t2_folder):
    "get data from all tier2 files from same run within a directory"

    
    run1_files = glob.glob(t2_folder + "/*run0001*.h5")

    file_list = []
    for filename in run1_files:
        df = pd.read_hdf(filename, "data")
        file_list.append(df)
    
    df_total = pd.concat(file_list, axis=0, ignore_index=True)
    keys = df_total.keys()
    data = df_total.to_numpy()

    #print(df_total.describe())
    #print(df_total.head())

    plt.figure()
    df_total.hist('e_ftp',  bins = 10000) #testing
    plt.ylabel("Frequency")
    plt.yscale("log")
    plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/plots/e_ftp_pandas.png")

    return keys, data


def obtain_key_data(data, keys, key, no_events):
    "obtain data values for a particular key, e.g. energy"

    key_index = (np.where(keys == key))[0][0]
    key_data = data[:no_events, key_index]

    return key_data

if __name__ =="__main__":
    main()